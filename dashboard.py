#!/usr/bin/env python3
"""Genomic ETL Dashboard with simplified dedup and filtering stats"""
from flask import Flask, render_template_string
import yaml
import sys
import boto3
from botocore.exceptions import BotoCoreError, ClientError
from botocore.config import Config

sys.path.insert(0, 'src')
from src.utils.database import Database

app = Flask(__name__)

with open('etl_config.yaml') as f:
    config = yaml.safe_load(f)

db = Database(config)

def _sum_s3_prefix(bucket, prefix, timeout=5):
    """Sum S3 prefix with timeout"""
    try:
        s3_config = Config(
            connect_timeout=timeout,
            read_timeout=timeout,
            retries={'max_attempts': 2}
        )
        s3 = boto3.client("s3", region_name=config["aws"].get("region", "us-east-1"), config=s3_config)
        paginator = s3.get_paginator("list_objects_v2")
        
        total_bytes = 0
        total_count = 0
        page_count = 0
        max_pages = 100  # Limit pages to prevent hanging
        
        for page in paginator.paginate(Bucket=bucket, Prefix=prefix, PaginationConfig={'MaxItems': 10000}):
            page_count += 1
            if page_count > max_pages:
                break
            for obj in page.get("Contents", []):
                total_bytes += int(obj.get("Size", 0))
                total_count += 1
        
        return total_bytes, total_count
    except Exception as e:
        print(f"S3 error for {prefix}: {e}")
        return None, None

def get_s3_storage_summary():
    """Get S3 storage summary with error handling"""
    try:
        bucket = config["aws"]["s3"]["bucket_name"]
        final_prefix = config["aws"]["s3"]["final_prefix"].rstrip("/") + "/"
        proteins_prefix = config["aws"]["s3"]["proteins_prefix"].rstrip("/") + "/"
        sources = config["sources"]["order"]
        
        per_source = []
        total_bytes = 0
        total_objects = 0
        
        for src in sources:
            final_src_prefix = f"{final_prefix}{src}/"
            proteins_src_prefix = f"{proteins_prefix}{src}/"
            
            b1, n1 = _sum_s3_prefix(bucket, final_src_prefix)
            b2, n2 = _sum_s3_prefix(bucket, proteins_src_prefix)
            
            if b1 is None or b2 is None:
                per_source.append((src, "Error", "Error"))
                continue
            
            src_bytes = b1 + b2
            src_objs = (n1 or 0) + (n2 or 0)
            total_bytes += src_bytes
            total_objects += src_objs
            per_source.append((src, src_objs, src_bytes))
        
        return total_bytes, per_source, total_objects
    except Exception as e:
        print(f"S3 summary error: {e}")
        return 0, [], 0

def _fmt_gb(b):
    if b is None or b == "Error":
        return "N/A"
    try:
        return f"{b/1_000_000_000:.2f} GB"
    except:
        return "N/A"

def _fmt_int(n):
    if n is None or n == "Error":
        return "N/A"
    try:
        return f"{n:,}"
    except:
        return "N/A"

DASHBOARD_HTML = """
<!DOCTYPE html>
<html>
<head>
    <title>Genomic ETL Dashboard</title>
    <meta http-equiv="refresh" content="30">
    <style>
        body { font-family: Arial; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
               padding: 20px; color: #333; }
        .container { max-width: 1300px; margin: 0 auto; }
        h1 { color: white; text-align: center; }
        .card { background: white; padding: 20px; margin: 20px 0; border-radius: 10px;
                box-shadow: 0 5px 15px rgba(0,0,0,0.2); }
        table { width: 100%; border-collapse: collapse; }
        th, td { padding: 10px 12px; text-align: left; border-bottom: 1px solid #eee; }
        th { background: #f7f7f7; color: #4b5bdc; font-weight: 600; }
        .stat { display: inline-block; margin: 10px 20px; }
        .stat-label { color: #666; font-size: 0.9em; }
        .stat-value { font-size: 1.8em; font-weight: bold; color: #667eea; }
        .flex { display: flex; gap: 20px; flex-wrap: wrap; }
        .col { flex: 1; min-width: 300px; }
        .muted { color: #777; font-size: 0.9em; }
        .kicker { font-size: 0.95em; color: #666; margin-top: 6px; }
        .center { text-align: center; color: white; margin-top: 20px; }
        .nowrap { white-space: nowrap; }
        .stack { display: flex; gap: 20px; flex-wrap: wrap; }
        .pill { display:inline-block; padding:6px 10px; background:#eef; border-radius:999px; color:#334; font-weight:600; }
        .success { color: #10b981; }
        .warning { color: #f59e0b; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Genomic ETL Pipeline Dashboard</h1>

        <div class="card">
            <h2>📊 Summary Statistics</h2>
            <div class="stat">
                <div class="stat-label">Total Entries</div>
                <div class="stat-value">{{ summary_stats.total_entries }}</div>
            </div>
            <div class="stat">
                <div class="stat-label">Total Base Pairs</div>
                <div class="stat-value">{{ total_bp }}</div>
            </div>
            <div class="stat">
                <div class="stat-label">Unique Species</div>
                <div class="stat-value">{{ total_species }}</div>
            </div>
            <div class="kicker">Includes both genomes and metagenomes; base pairs shown in bp.</div>
        </div>

        <div class="card">
            <h2>🔄 Deduplication Statistics</h2>
            <div class="stack">
                <div class="pill">Total Entries: {{ "{:,}".format(dedup_stats.total_entries) }}</div>
                <div class="pill success">Unique Entries: {{ "{:,}".format(dedup_stats.unique_entries) }}</div>
                <div class="pill warning">Duplicates: {{ "{:,}".format(dedup_stats.duplicate_entries) }}</div>
            </div>
            <div class="kicker">
                Uniqueness rate: {{ "%.2f"|format((dedup_stats.unique_entries / dedup_stats.total_entries * 100) if dedup_stats.total_entries > 0 else 0) }}%
                • Duplicates are identified by sequence hash or accession match
            </div>
        </div>

        <div class="card">
            <h2>🔍 Filtering Statistics</h2>
            <div class="stack">
                <div class="pill">Total Contigs: {{ "{:,}".format(filtering_stats.total_contigs) }}</div>
                <div class="pill success">Contigs Kept: {{ "{:,}".format(filtering_stats.contigs_kept) }}</div>
                <div class="pill warning">Contigs Removed: {{ "{:,}".format(filtering_stats.contigs_removed) }}</div>
            </div>
            <div class="kicker">
                Retention rate: {{ "%.2f"|format((filtering_stats.contigs_kept / filtering_stats.total_contigs * 100) if filtering_stats.total_contigs > 0 else 0) }}%
                • Metagenome contigs < {{ min_contig_length }} bp are filtered out
            </div>
        </div>

        <div class="card">
            <h2>🗄️ Database Statistics (by Source)</h2>
            <table>
                <thead>
                    <tr>
                        <th>Source</th>
                        <th>Entries</th>
                        <th>Base Pairs (bp)</th>
                        <th>Species</th>
                    </tr>
                </thead>
                <tbody>
                    {% for stat in stats %}
                    <tr>
                        <td>{{ stat[0] }}</td>
                        <td>{{ "{:,}".format(stat[1]) }}</td>
                        <td>{{ "{:,}".format(stat[2]) }}</td>
                        <td>{{ stat[3] }}</td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
        </div>

        <div class="flex">
            <div class="card col">
                <h2>🧭 Counts by Kingdom (Genomes)</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Kingdom</th>
                            <th>Entries</th>
                            <th>Base Pairs (bp)</th>
                            <th>Species</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for k in by_kingdom %}
                        <tr>
                            <td>{{ k[0] }}</td>
                            <td>{{ "{:,}".format(k[1]) }}</td>
                            <td>{{ "{:,}".format(k[2]) }}</td>
                            <td>{{ k[3] }}</td>
                        </tr>
                        {% endfor %}
                        {% if not by_kingdom %}
                        <tr><td colspan="4" class="muted">No kingdom data yet.</td></tr>
                        {% endif %}
                    </tbody>
                </table>
            </div>

            <div class="card col">
                <h2>🌍 Counts by Origin (Metagenomes)</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Origin</th>
                            <th>Entries</th>
                            <th>Base Pairs (bp)</th>
                            <th>Species</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for o in by_origin %}
                        <tr>
                            <td>{{ o[0] }}</td>
                            <td>{{ "{:,}".format(o[1]) }}</td>
                            <td>{{ "{:,}".format(o[2]) }}</td>
                            <td>{{ o[3] }}</td>
                        </tr>
                        {% endfor %}
                        {% if not by_origin %}
                        <tr><td colspan="4" class="muted">No origin data yet.</td></tr>
                        {% endif %}
                    </tbody>
                </table>
            </div>
        </div>

        <div class="card">
            <h2>💾 S3 Storage</h2>
            <div class="stack">
                <div class="pill">Total: {{ total_s3_gb }}</div>
                <div class="pill">Objects: {{ total_s3_objects }}</div>
            </div>
            <table style="margin-top:12px;">
                <thead>
                    <tr>
                        <th>Source</th>
                        <th>Objects</th>
                        <th>Size (GB)</th>
                    </tr>
                </thead>
                <tbody>
                    {% for row in s3_by_source %}
                    <tr>
                        <td>{{ row[0] }}</td>
                        <td>{{ row[1] }}</td>
                        <td>{{ row[2] }}</td>
                    </tr>
                    {% endfor %}
                    {% if not s3_by_source %}
                    <tr><td colspan="3" class="muted">No S3 data available</td></tr>
                    {% endif %}
                </tbody>
            </table>
            <div class="muted">Totals include both genome (final/) and protein (proteins/) prefixes. S3 stats may be approximate due to pagination limits.</div>
        </div>

        <div class="center">Auto-refresh every 30 seconds | Last updated: {{ timestamp }}</div>
    </div>
</body>
</html>
"""

@app.route('/')
def index():
    from datetime import datetime
    
    try:
        stats = db.get_stats()
        by_kingdom = db.get_counts_by_kingdom()
        by_origin = db.get_counts_by_origin()
        dedup_stats = db.get_dedup_stats()
        filtering_stats = db.get_filtering_stats()

        total_entries = sum(s[1] for s in stats)
        total_bp = sum(s[2] for s in stats)
        total_species = sum(s[3] for s in stats)

        # S3 stats with timeout
        total_bytes, per_source, total_objects = get_s3_storage_summary()
        s3_by_source = [(src, _fmt_int(objs), _fmt_gb(bytes_)) for (src, objs, bytes_) in per_source]
        total_s3_gb = _fmt_gb(total_bytes)
        total_s3_objects = _fmt_int(total_objects)

        min_contig_length = config.get('filtering', {}).get('metagenomes', {}).get('min_contig_length', 2000)

        return render_template_string(
            DASHBOARD_HTML,
            stats=stats,
            by_kingdom=by_kingdom,
            by_origin=by_origin,
            dedup_stats=dedup_stats,
            filtering_stats=filtering_stats,
            summary_stats={'total_entries': f"{total_entries:,}"},
            total_bp=f"{total_bp:,} bp",
            total_species=total_species,
            total_s3_gb=total_s3_gb,
            total_s3_objects=total_s3_objects,
            s3_by_source=s3_by_source,
            min_contig_length=min_contig_length,
            timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        )
    except Exception as e:
        return f"<h1>Dashboard Error</h1><pre>{str(e)}</pre>", 500

if __name__ == '__main__':
    port = config['dashboard']['port']
    print(f"🚀 Dashboard at http://localhost:{port}")
    print(f"   Access from browser: http://<your-ec2-ip>:{port}")
    app.run(host='0.0.0.0', port=port, debug=False, threaded=True)