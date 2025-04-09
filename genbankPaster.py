import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict

# Gene function categorization dictionary
GENE_FUNCTIONS = {
    # Metabolic genes
    'metabolism': ['dehydrogenase', 'synthase', 'kinase', 'transferase', 'reductase', 
                  'oxidase', 'oxygenase', 'carboxylase', 'decarboxylase', 'isomerase',
                  'phosphatase', 'hydrolase', 'polymerase', 'ligase', 'lyase'],
    
    # Structural genes
    'structural': ['actin', 'tubulin', 'collagen', 'keratin', 'myosin', 'fibrin',
                  'laminin', 'fibronectin', 'elastin', 'spectrin', 'dystrophin'],
    
    # Transport genes
    'transport': ['transporter', 'channel', 'porin', 'pump', 'carrier', 'symporter',
                 'antiporter', 'permease', 'import', 'export', 'ATPase'],
    
    # Regulatory genes
    'regulatory': ['transcription factor', 'repressor', 'activator', 'enhancer',
                  'silencer', 'operator', 'promoter', 'regulator', 'homeobox'],
    
    # Signaling genes
    'signaling': ['receptor', 'G protein', 'kinase', 'phosphatase', 'adapter',
                 'second messenger', 'calmodulin', 'ras', 'mapk', 'jak', 'stat'],
    
    # Immune genes
    'immune': ['antibody', 'immunoglobulin', 'MHC', 'histocompatibility', 'interleukin',
              'interferon', 'complement', 'defensin', 'toll-like', 'lysozyme'],
    
    # Developmental genes
    'developmental': ['homeobox', 'notch', 'wingless', 'hedgehog', 'bone morphogenetic',
                     'fibroblast growth', 'transforming growth', 'sonic hedgehog'],
    
    # Housekeeping genes
    'housekeeping': ['GAPDH', 'actin', 'tubulin', 'ubiquitin', 'ribosomal protein',
                    'histone', 'elongation factor', 'RNA polymerase'],
    
    # Unknown/other
    'unknown': ['unknown', 'hypothetical', 'uncharacterized', 'putative']
}

def categorize_gene(gene_name, product_description):
    """Categorize gene based on its name and product description."""
    search_text = f"{gene_name.lower()} {product_description.lower()}"
    
    for category, keywords in GENE_FUNCTIONS.items():
        if any(keyword.lower() in search_text for keyword in keywords):
            return category
    return 'other'

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
                
                # Categorize the gene
                category = categorize_gene(gene_name, product)
                
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
                    species,
                    category
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
    gene_categories = defaultdict(int)
    
    if os.path.isdir(input_path):
        files = [f for f in os.listdir(input_path) 
                if f.lower().endswith((".gb", ".genbank", ".gbk"))]
    else:
        files = [input_path]
    
    for filepath in files:
        if os.path.isdir(input_path):
            filepath = os.path.join(input_path, filepath)
        
        features = extract_feature_from_genbank(filepath, feature_type, feature_names)
        for gene_name, seq_record, species, category in features:
            # Skip exact duplicates
            if seq_record.id in seen_ids:
                continue
                
            extracted_data[gene_name].append((seq_record, category))
            all_genes.add(gene_name)
            all_species.add(species)
            seen_ids.add(seq_record.id)
            gene_categories[category] += 1
    
    return extracted_data, all_genes, all_species, gene_categories

def export_to_fasta(output_dir, extracted_data):
    """Exports sequences to FASTA files with improved naming."""
    os.makedirs(output_dir, exist_ok=True)
    for gene_name, seq_records in extracted_data.items():
        clean_name = "".join(c if c.isalnum() else "_" for c in gene_name)
        fasta_path = os.path.join(output_dir, f"{clean_name}.fasta")
        with open(fasta_path, "w") as f:
            SeqIO.write([sr[0] for sr in seq_records], f, "fasta")
        print(f"Exported {len(seq_records)} sequences to {fasta_path}")

def generate_html_report(output_dir, extracted_data, all_genes, all_species, gene_categories):
    """Generates comprehensive HTML report with species information and phylogenetic gene suggestions."""
    # Statistics calculation
    gene_stats = {gene: len(seqs) for gene, seqs in extracted_data.items()}
    total_sequences = sum(gene_stats.values())
    sorted_genes = sorted(gene_stats.items(), key=lambda x: (-x[1], x[0]))
    
    # Prepare category data
    category_data = []
    for category, count in gene_categories.items():
        category_data.append({
            'name': category.capitalize(),
            'count': count,
            'percentage': count/total_sequences,
            'genes': defaultdict(int)
        })
    
    # Count genes per category
    for gene, seq_records in extracted_data.items():
        category = seq_records[0][1]  # Get category from first record
        for cat in category_data:
            if cat['name'].lower() == category:
                cat['genes'][gene] += len(seq_records)
                break
    
    # Sort categories by count
    category_data.sort(key=lambda x: -x['count'])
    
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
        f"<tr><td>{sp}</td><td>{sum(1 for seqs in extracted_data.values() for s in seqs if sp in s[0].description)}</td></tr>"
        for sp in sorted(all_species)
    )
    
    sequences_table = "\n".join(
        f"""<tr>
            <td>{gene}</td>
            <td>{seq.id}</td>
            <td>{seq.description.split('|')[-1].strip()}</td>
            <td>{len(seq.seq)}</td>
            <td>{category}</td>
        </tr>"""
        for gene, seq_records in extracted_data.items()
        for seq, category in seq_records
    )
    
    category_chart_data = ",".join(
        f"{{label: '{cat['name']}', value: {cat['count']}}}"
        for cat in category_data
    )
    
    category_table = "\n".join(
        f"""<tr>
            <td>{cat['name']}</td>
            <td>{cat['count']}</td>
            <td>{cat['percentage']:.1%}</td>
            <td>{', '.join(sorted(cat['genes'].keys())[:5]) + ('...' if len(cat['genes']) > 5 else '')}</td>
        </tr>"""
        for cat in category_data
    )
    
    # Phylogenetic marker gene suggestions
    PHYLOGENETIC_MARKERS = {
        'plastome': {
            'Highly Recommended': ['rbcL', 'matK', 'ndhF', 'atpB', 'rpoC2'],
            'Recommended': ['psbA', 'psbB', 'psbC', 'psbD', 'rpoB', 'ycf1', 'ycf2'],
            'rRNA genes': ['rrn16', 'rrn23', 'rrn5', 'rrn4.5'],
            'tRNA genes': ['trnK-UUU', 'trnL-UAA', 'trnF-GAA', 'trnH-GUG']
        },
        'mitochondrial': {
            'Protein-coding': ['cox1', 'cox2', 'cox3', 'cob', 'nad1', 'nad2', 'nad3', 'nad4', 'nad5', 'nad6', 'atp6', 'atp8'],
            'rRNA genes': ['rrnL', 'rrnS'],
            'tRNA genes': ['trnM', 'trnW', 'trnQ', 'trnY']
        },
        'bacterial': {
            'Universal': ['16S rRNA', '23S rRNA'],
            'Protein-coding': ['rpoB', 'gyrB', 'recA', 'dnaK', 'tuf', 'fusA'],
            'Housekeeping': ['rpoD', 'groEL', 'atpD', 'dnaJ']
        }
    }
    
    # Prepare phylogenetic markers section
    marker_html = ""
    for genome_type, categories in PHYLOGENETIC_MARKERS.items():
        marker_html += f"""
        <div class="marker-category">
            <h3>{genome_type.capitalize()} Markers</h3>
            <table>
                <tr><th>Category</th><th>Suggested Genes</th><th>In Your Data</th></tr>
        """
        
        for category, genes in categories.items():
            present_genes = []
            missing_genes = []
            
            for gene in genes:
                if gene.lower() in (g.lower() for g in all_genes):
                    present_genes.append(gene)
                else:
                    missing_genes.append(gene)
            
            present_html = "<span style='color:green'>" + ", ".join(present_genes) + "</span>" if present_genes else "-"
            missing_html = "<span style='color:red'>" + ", ".join(missing_genes) + "</span>" if missing_genes else "-"
            
            all_genes_list = present_genes + missing_genes
            marker_html += f"""
                <tr>
                    <td>{category}</td>
                    <td>{", ".join(all_genes_list)}</td>
                    <td>{present_html if present_genes else "None"} 
                        {f"<br>(Missing: {missing_html})" if missing_genes else ""}</td>
                </tr>
            """
        
        marker_html += """
            </table>
        </div>
        """
    
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
        #categoryChart {{ width: 100%; height: 400px; margin: 20px 0; }}
        .two-column {{ display: flex; gap: 20px; }}
        .column {{ flex: 1; }}
        .marker-category {{
            margin-bottom: 20px;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 5px;
        }}
        .marker-category h3 {{
            margin-top: 0;
            color: #2E86C1;
        }}
        .advice {{
            background: #e8f4f8;
            padding: 15px;
            border-radius: 5px;
            margin-top: 20px;
        }}
    </style>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
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
        <h2>Functional Categories</h2>
        <div class="two-column">
            <div class="column">
                <canvas id="categoryChart"></canvas>
            </div>
            <div class="column">
                <table>
                    <tr><th>Category</th><th>Count</th><th>Percentage</th><th>Example Genes</th></tr>
                    {category_table}
                </table>
            </div>
        </div>
    </div>
    
    <div class="section">
        <h2>Gene Statistics</h2>
        <div class="tabs">
            <div class="tab active" onclick="showTab('common-genes')">Most Common</div>
            <div class="tab" onclick="showTab('all-genes')">All Genes</div>
            <div class="tab" onclick="showTab('phylogenetic-markers')">Phylogenetic Markers</div>
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
        
        <div id="phylogenetic-markers" class="tab-content">
            <h3>Suggested Genes for Phylogenetic Analysis</h3>
            <p>Based on your extracted genes and common phylogenetic markers, here are recommendations 
            for genes to use in your phylogenetic analysis:</p>
            
            {marker_html}
            
            <div class="advice">
                <h4>Phylogenetic Analysis Advice:</h4>
                <ul>
                    <li><strong>For plastomes:</strong> The combination of <em>rbcL + matK</em> is commonly used as a DNA barcode for plants. For deeper phylogeny, consider adding <em>ndhF</em> or <em>rpoC2</em>.</li>
                    <li><strong>For mitochondrial genomes:</strong> <em>cox1</em> is the standard animal barcode. For deeper analysis, combine with <em>cob</em> and <em>nad</em> genes.</li>
                    <li><strong>For bacteria:</strong> <em>16S rRNA</em> is standard for broad classification. For finer resolution, use protein-coding genes like <em>rpoB</em> or <em>gyrB</em>.</li>
                    <li>When possible, use concatenated datasets of multiple genes for better resolution.</li>
                    <li>For organelle genomes, consider using whole genome alignments if you have enough conserved regions.</li>
                    <li>For closely related species, non-coding regions may provide better resolution.</li>
                    <li>Always verify that your selected markers have appropriate evolutionary rates for your timescale.</li>
                </ul>
            </div>
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
                <th>Category</th>
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
        
        // Draw category chart
        document.addEventListener('DOMContentLoaded', function() {{
            const ctx = document.getElementById('categoryChart').getContext('2d');
            const data = [{category_chart_data}];
            
            const backgroundColors = [
                'rgba(255, 99, 132, 0.7)',
                'rgba(54, 162, 235, 0.7)',
                'rgba(255, 206, 86, 0.7)',
                'rgba(75, 192, 192, 0.7)',
                'rgba(153, 102, 255, 0.7)',
                'rgba(255, 159, 64, 0.7)',
                'rgba(199, 199, 199, 0.7)',
                'rgba(83, 102, 255, 0.7)',
                'rgba(255, 99, 255, 0.7)'
            ];
            
            new Chart(ctx, {{
                type: 'doughnut',
                data: {{
                    labels: data.map(item => item.label),
                    datasets: [{{
                        data: data.map(item => item.value),
                        backgroundColor: backgroundColors,
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    plugins: {{
                        legend: {{
                            position: 'right',
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const label = context.label || '';
                                    const value = context.raw || 0;
                                    const total = context.dataset.data.reduce((a, b) => a + b, 0);
                                    const percentage = Math.round((value / total) * 100);
                                    return `${{label}}: ${{value}} (${{percentage}}%)`;
                                }}
                            }}
                        }}
                    }}
                }}
            }});
        }});
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
    
    extracted_data, all_genes, all_species, gene_categories = process_files(
        args.input_path,
        args.feature_type,
        [name.lower() for name in args.feature_names] if args.feature_names else None
    )
    
    os.makedirs(args.output_dir, exist_ok=True)
    export_to_fasta(args.output_dir, extracted_data)
    generate_html_report(args.output_dir, extracted_data, all_genes, all_species, gene_categories)

if __name__ == "__main__":
    main()
