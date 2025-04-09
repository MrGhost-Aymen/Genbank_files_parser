import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict

def sanitize_id(original_id):
    """Sanitizes sequence IDs for compatibility with phylogenetic software."""
    illegal_chars = [" ", ":", ",", ")", "(", ";", "]", "[", "'", '"', "|", "."]
    sanitized_id = original_id
    for char in illegal_chars:
        sanitized_id = sanitized_id.replace(char, "_")
    sanitized_id = "".join([c if c.isalnum() or c == "_" else "" for c in sanitized_id])
    return sanitized_id or "Unknown"

def extract_feature_from_genbank(genbank_file, feature_type, feature_names=None):
    """Extracts features from GenBank files with enhanced metadata handling."""
    extracted_features = []
    try:
        record = SeqIO.read(genbank_file, "genbank")
        accession = record.id
        species = record.description.split("[")[0].strip() if "[" in record.description else record.description
        
        for feature in record.features:
            if feature.type == feature_type:
                gene_names = feature.qualifiers.get("gene", ["Unknown"])
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                
                if feature_names:  # Check if we should extract this feature
                    should_extract = False
                    for name in feature_names:
                        # Check against both gene names and product names
                        if (any(name.lower() in gn.lower() for gn in gene_names) or
                            any(name.lower() in pn.lower() for pn in feature.qualifiers.get("product", [""]))):
                            should_extract = True
                            break
                    if not should_extract:
                        continue
                
                # Get the primary gene name
                gene_name = gene_names[0] if gene_names[0] != "Unknown" else product.split()[0]
                
                # Extract sequence and create record
                seq = feature.location.extract(record.seq)
                seq_id = f"{accession}_{gene_name}"
                seq_id = sanitize_id(seq_id)
                
                extracted_features.append((
                    gene_name,
                    SeqRecord(
                        Seq(str(seq)),
                        id=seq_id,
                        description=f"{product} | {species}"
                    ),
                    species
                ))
    except Exception as e:
        print(f"Error processing {genbank_file}: {e}")
    return extracted_features

def process_files(input_path, feature_type, feature_names=None):
    """Processes GenBank files with duplicate gene handling."""
    extracted_data = defaultdict(list)
    all_genes = set()
    all_species = set()
    seen_ids = set()
    
    if os.path.isdir(input_path):
        files = [f for f in os.listdir(input_path) 
                if f.lower().endswith((".gb", ".genbank", ".gbk"))]
    else:
        files = [input_path]
    
    for filepath in files:
        if os.path.isdir(input_path):
            filepath = os.path.join(input_path, filepath)
        
        features = extract_feature_from_genbank(filepath, feature_type, feature_names)
        for gene_name, seq_record, species in features:
            # Skip exact duplicates
            if seq_record.id in seen_ids:
                continue
                
            extracted_data[gene_name].append(seq_record)
            all_genes.add(gene_name)
            all_species.add(species)
            seen_ids.add(seq_record.id)
    
    return extracted_data, all_genes, all_species

def export_to_fasta(output_dir, extracted_data):
    """Exports sequences to FASTA files with improved naming."""
    os.makedirs(output_dir, exist_ok=True)
    for gene_name, seq_records in extracted_data.items():
        clean_name = "".join(c if c.isalnum() else "_" for c in gene_name)
        fasta_path = os.path.join(output_dir, f"{clean_name}.fasta")
        with open(fasta_path, "w") as f:
            SeqIO.write(seq_records, f, "fasta")
        print(f"Exported {len(seq_records)} sequences to {fasta_path}")

def generate_html_report(output_dir, extracted_data, all_genes, all_species):
    """Generates comprehensive HTML report with species information."""
    # Statistics calculation
    gene_stats = {gene: len(seqs) for gene, seqs in extracted_data.items()}
    total_sequences = sum(gene_stats.values())
    sorted_genes = sorted(gene_stats.items(), key=lambda x: (-x[1], x[0]))
    
    # Prepare report sections
    common_genes = "\n".join(
        f"<tr><td>{gene}</td><td>{count}</td><td>{count/total_sequences:.1%}</td></tr>"
        for gene, count in sorted_genes[:20]
    )
    
    all_genes_list = "\n".join(
        f"<tr><td>{gene}</td><td>{count}</td><td>{count/total_sequences:.1%}</td></tr>"
        for gene, count in sorted_genes
    )
    
    species_list = "\n".join(
        f"<tr><td>{sp}</td><td>{sum(1 for seqs in extracted_data.values() for s in seqs if sp in s.description)}</td></tr>"
        for sp in sorted(all_species)
    )
    
    sequences_table = "\n".join(
        f"""<tr>
            <td>{gene}</td>
            <td>{seq.id}</td>
            <td>{seq.description.split('|')[-1].strip()}</td>
            <td>{len(seq.seq)}</td>
        </tr>"""
        for gene, seqs in extracted_data.items()
        for seq in seqs
    )
    
    html_content = f"""
<html>
<head>
    <title>GenBank Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1, h2 {{ color: #2E86C1; }}
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; position: sticky; top: 0; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
        .section {{ margin-bottom: 30px; padding: 15px; border-radius: 5px; 
                  box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .stats {{ display: flex; gap: 20px; margin-bottom: 20px; }}
        .stat-box {{ flex: 1; padding: 15px; background: #f8f9fa; 
                    border-radius: 5px; text-align: center; }}
        .stat-value {{ font-size: 24px; font-weight: bold; color: #2E86C1; }}
        .tabs {{ display: flex; gap: 5px; margin-bottom: 10px; }}
        .tab {{ padding: 10px 15px; background: #e9ecef; cursor: pointer; 
               border-radius: 5px 5px 0 0; }}
        .tab.active {{ background: #2E86C1; color: white; }}
        .tab-content {{ display: none; }}
        .tab-content.active {{ display: block; }}
    </style>
</head>
<body>
    <h1>GenBank Feature Extraction Report</h1>
    
    <div class="stats">
        <div class="stat-box">
            <div class="stat-value">{len(all_genes)}</div>
            <div>Unique Genes</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{total_sequences}</div>
            <div>Total Sequences</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{len(all_species)}</div>
            <div>Species</div>
        </div>
    </div>
    
    <div class="section">
        <h2>Gene Statistics</h2>
        <div class="tabs">
            <div class="tab active" onclick="showTab('common-genes')">Most Common</div>
            <div class="tab" onclick="showTab('all-genes')">All Genes</div>
        </div>
        
        <div id="common-genes" class="tab-content active">
            <table>
                <tr><th>Gene</th><th>Count</th><th>Percentage</th></tr>
                {common_genes}
            </table>
        </div>
        
        <div id="all-genes" class="tab-content">
            <table>
                <tr><th>Gene</th><th>Count</th><th>Percentage</th></tr>
                {all_genes_list}
            </table>
        </div>
    </div>
    
    <div class="section">
        <h2>Species Information</h2>
        <table>
            <tr><th>Species</th><th>Sequence Count</th></tr>
            {species_list}
        </table>
    </div>
    
    <div class="section">
        <h2>Extracted Sequences</h2>
        <table>
            <tr>
                <th>Gene</th>
                <th>Accession</th>
                <th>Species</th>
                <th>Length</th>
            </tr>
            {sequences_table}
        </table>
    </div>
    
    <script>
        function showTab(tabId) {{
            document.querySelectorAll('.tab-content').forEach(tab => {{
                tab.classList.remove('active');
            }});
            document.querySelectorAll('.tab').forEach(tab => {{
                tab.classList.remove('active');
            }});
            document.getElementById(tabId).classList.add('active');
            event.currentTarget.classList.add('active');
        }}
    </script>
</body>
</html>
    """
    
    report_path = os.path.join(output_dir, "detailed_report.html")
    with open(report_path, "w") as f:
        f.write(html_content)
    print(f"Generated comprehensive report at {report_path}")

def main():
    parser = argparse.ArgumentParser(description="Advanced GenBank feature extractor")
    parser.add_argument("input_path", help="GenBank file or directory")
    parser.add_argument("--feature_type", required=True, 
                       choices=["CDS", "tRNA", "rRNA", "gene", "mRNA"],
                       help="Type of feature to extract")
    parser.add_argument("--feature_names", nargs="+", 
                       help="Specific gene/product names to extract")
    parser.add_argument("--output_dir", default="genbank_output",
                       help="Output directory (default: genbank_output)")
    
    args = parser.parse_args()
    
    extracted_data, all_genes, all_species = process_files(
        args.input_path,
        args.feature_type,
        [name.lower() for name in args.feature_names] if args.feature_names else None
    )
    
    os.makedirs(args.output_dir, exist_ok=True)
    export_to_fasta(args.output_dir, extracted_data)
    generate_html_report(args.output_dir, extracted_data, all_genes, all_species)

if __name__ == "__main__":
    main()
