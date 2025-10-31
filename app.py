import streamlit as st
import json
import requests
from typing import Dict, List, Tuple
import time
from io import StringIO
import pandas as pd
import plotly.express as px

# Ensembl REST API base URL
ENSEMBL_API_BASE = "https://rest.ensembl.org"

st.set_page_config(
    page_title="Variant to Gene Annotation", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        margin-bottom: 0.5rem;
    }
    .section-header {
        font-size: 1.5rem;
        font-weight: bold;
        color: #2c3e50;
        margin-top: 2rem;
        margin-bottom: 1rem;
        border-bottom: 2px solid #3498db;
        padding-bottom: 0.5rem;
    }
    .metric-container {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #3498db;
    }
    .success-box {
        background-color: #d4edda;
        border: 1px solid #c3e6cb;
        border-radius: 0.5rem;
        padding: 1rem;
        color: #155724;
    }
    .warning-box {
        background-color: #fff3cd;
        border: 1px solid #ffc107;
        border-radius: 0.5rem;
        padding: 1rem;
        color: #856404;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<h1 class="main-header">Variant to Gene Annotation Tool</h1>', unsafe_allow_html=True)
st.markdown("Analyze SNPs and annotate with nearest genes using Ensembl REST API")

# Load and parse the CSV file
@st.cache_data
def load_variants(file_path: str) -> pd.DataFrame:
    """Load and parse the variants CSV file."""
    try:
        df = pd.read_csv(file_path)
        return df
    except Exception as e:
        st.error(f"Error loading CSV file: {e}")
        return pd.DataFrame()

def validate_variants(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Identify duplicated and malformed entries."""
    # Identify duplicates
    duplicates = df[df.duplicated(subset=['variant_id'], keep=False)]
    
    # Identify malformed entries
    malformed_list = []
    
    # Check for missing values in required columns
    required_cols = ['variant_id', 'chrom', 'pos', 'ref', 'alt']
    missing_data = df[df[required_cols].isnull().any(axis=1)]
    if not missing_data.empty:
        malformed_list.append(missing_data)
    
    # Check for invalid chromosome values (should be 1-22, X, Y, MT)
    valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    invalid_chrom = df[~df['chrom'].astype(str).isin(valid_chroms)]
    if not invalid_chrom.empty:
        malformed_list.append(invalid_chrom)
    
    # Check for invalid position values (should be positive integers)
    invalid_pos = df[(df['pos'] <= 0) | (~df['pos'].apply(lambda x: isinstance(x, (int, float)) and x > 0))]
    if not invalid_pos.empty:
        malformed_list.append(invalid_pos)
    
    # Check for invalid ref/alt values (should be valid nucleotide codes)
    valid_bases = ['A', 'T', 'G', 'C']
    invalid_ref = df[~df['ref'].isin(valid_bases)]
    if not invalid_ref.empty:
        malformed_list.append(invalid_ref)
    
    invalid_alt = df[~df['alt'].isin(valid_bases)]
    if not invalid_alt.empty:
        malformed_list.append(invalid_alt)
    
    # Combine all malformed entries and remove duplicates
    if malformed_list:
        malformed = pd.concat(malformed_list).drop_duplicates()
    else:
        malformed = pd.DataFrame()
    
    return duplicates, malformed

def get_nearest_gene(chrom: str, pos: int, window: int = 50000) -> Dict:
    """Get the nearest gene for a variant using Ensembl REST API."""
    # Convert chromosome format (remove 'chr' prefix if present)
    chrom_clean = str(chrom).replace('chr', '')
    
    # Define the region window around the variant
    start = max(1, pos - window)
    end = pos + window
    
    # Build the API endpoint
    endpoint = f"{ENSEMBL_API_BASE}/overlap/region/human/{chrom_clean}:{start}-{end}"
    params = {
        "feature": "gene",
        "content-type": "application/json"
    }
    
    try:
        response = requests.get(endpoint, params=params, headers={"Content-Type": "application/json"})
        response.raise_for_status()
        
        genes = response.json()
        
        if not genes:
            return {"gene_id": None, "gene_name": None, "distance": None, "strand": None}
        
        # Find the nearest gene
        nearest_gene = None
        min_distance = float('inf')
        
        for gene in genes:
            gene_start = gene.get('start', 0)
            gene_end = gene.get('end', 0)
            
            # Calculate distance to gene
            if pos < gene_start:
                distance = gene_start - pos
            elif pos > gene_end:
                distance = pos - gene_end
            else:
                distance = 0  # Variant is within the gene
            
            if distance < min_distance:
                min_distance = distance
                nearest_gene = gene
        
        if nearest_gene:
            return {
                "gene_id": nearest_gene.get('id', 'N/A'),
                "gene_name": nearest_gene.get('external_name', nearest_gene.get('id', 'N/A')),
                "distance": min_distance,
                "strand": nearest_gene.get('strand', None)
            }
        else:
            return {"gene_id": None, "gene_name": None, "distance": None, "strand": None}
    
    except requests.exceptions.RequestException as e:
        st.warning(f"API error for {chrom}:{pos}: {str(e)}")
        return {"gene_id": None, "gene_name": None, "distance": None, "strand": None}

def annotate_variants(df: pd.DataFrame, progress_bar) -> pd.DataFrame:
    """Annotate variants with nearest genes."""
    total_variants = len(df)
    
    # Collect annotations
    annotations = []
    
    # Iterate over rows
    for idx, row in df.iterrows():
        gene_info = get_nearest_gene(row['chrom'], int(row['pos']))
        
        annotations.append({
            'nearest_gene_id': gene_info['gene_id'],
            'nearest_gene_name': gene_info['gene_name'],
            'distance_to_gene': gene_info['distance'],
            'gene_strand': gene_info['strand']
        })
        
        # Update progress bar
        progress_bar.progress((idx + 1) / total_variants)
        
        # Be respectful to the API - add a small delay
        time.sleep(0.1)
    
    # Create annotation DataFrame and join with original
    annotation_df = pd.DataFrame(annotations)
    annotated_df = pd.concat([df.reset_index(drop=True), annotation_df], axis=1)
    
    return annotated_df

def generate_summary(annotated_df: pd.DataFrame) -> Dict:
    """Generate summary JSON with chromosome counts and top 10 genes."""
    # Counts per chromosome
    chrom_counts = annotated_df['chrom'].value_counts().to_dict()
    chrom_counts = {str(k): int(v) for k, v in chrom_counts.items()}
    
    # Top 10 genes with the most variants (excluding null values)
    gene_counts = annotated_df[annotated_df['nearest_gene_name'].notna()]['nearest_gene_name'].value_counts().head(10)
    top_10_genes = {str(k): int(v) for k, v in gene_counts.items()}
    
    # Count variants with/without annotations
    variants_with_annotation = annotated_df['nearest_gene_name'].notna().sum()
    variants_without_annotation = annotated_df['nearest_gene_name'].isna().sum()
    
    summary = {
        "chromosome_counts": chrom_counts,
        "top_10_genes": top_10_genes,
        "total_variants": len(annotated_df),
        "variants_with_gene_annotation": int(variants_with_annotation),
        "variants_without_gene_annotation": int(variants_without_annotation)
    }
    
    return summary

# Main app logic
file_path = "sample_variants.csv"

# Load data
df = load_variants(file_path)

if not df.empty:
    # Create tabs for better organization
    tab1, tab2, tab3, tab4 = st.tabs(["Overview", "Validation", "Annotation", "Summary"])
    
    with tab1:
        st.markdown('<h2 class="section-header">Data Overview</h2>', unsafe_allow_html=True)
        
        # Enhanced metrics with styling
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Variants", len(df), help="Total number of variants in the dataset")
        with col2:
            st.metric("Unique Variants", df['variant_id'].nunique(), help="Number of unique variant IDs")
        with col3:
            st.metric("Chromosomes", df['chrom'].nunique(), help="Number of different chromosomes")
        with col4:
            st.metric("Unique Positions", df[['chrom', 'pos']].drop_duplicates().shape[0], help="Number of unique chromosome-position combinations")
        
        # Display raw data with better formatting
        st.markdown('<h3 style="margin-top: 2rem;">Raw Data Preview</h3>', unsafe_allow_html=True)
        
        st.dataframe(
            df.head(20),
            width='stretch',
            hide_index=True,
            height=400
        )
        
        if len(df) > 20:
            st.caption(f"Showing first 20 of {len(df)} variants. Use annotation tab to see all variants with gene annotations.")
        
        # Chromosome distribution visualization
        st.markdown('<h3 style="margin-top: 2rem;">Chromosome Distribution</h3>', unsafe_allow_html=True)
        chrom_counts_df = df['chrom'].value_counts().sort_index().reset_index()
        chrom_counts_df.columns = ['chrom', 'count']
        
        fig_chrom = px.bar(
            chrom_counts_df,
            x='chrom',
            y='count',
            title='Variants per Chromosome',
            labels={'chrom': 'Chromosome', 'count': 'Number of Variants'},
            color='count',
            color_continuous_scale='Blues'
        )
        fig_chrom.update_layout(
            showlegend=False,
            height=400,
            xaxis_title="Chromosome",
            yaxis_title="Number of Variants"
        )
        st.plotly_chart(fig_chrom, width='stretch')
    
    with tab2:
        st.markdown('<h2 class="section-header">Data Validation</h2>', unsafe_allow_html=True)
        
        duplicates, malformed = validate_variants(df)
        
        # Duplicates section
        col1, col2 = st.columns([1, 3])
        with col1:
            if not duplicates.empty:
                st.error(f"**{len(duplicates)} duplicate entries found**")
            else:
                st.success("**No duplicates found**")
        
        if not duplicates.empty:
            with st.expander("View Duplicate Entries", expanded=True):
                st.dataframe(duplicates, width='stretch', hide_index=True)
        
        st.divider()
        
        # Malformed entries section
        col1, col2 = st.columns([1, 3])
        with col1:
            if not malformed.empty:
                st.warning(f"**{len(malformed)} malformed entries found**")
            else:
                st.success("**No malformed entries found**")
        
        if not malformed.empty:
            with st.expander("View Malformed Entries", expanded=True):
                st.dataframe(malformed, width='stretch', hide_index=True)
        
        # Data quality summary
        st.markdown('<h3 style="margin-top: 2rem;">Data Quality Summary</h3>', unsafe_allow_html=True)
        quality_metrics = {
            'Total Variants': len(df),
            'Duplicates': len(duplicates),
            'Malformed': len(malformed),
            'Valid Variants': len(df) - len(duplicates) - len(malformed),
            'Quality Score': f"{(len(df) - len(duplicates) - len(malformed)) / len(df) * 100:.1f}%"
        }
        
        quality_df = pd.DataFrame(list(quality_metrics.items()), columns=['Metric', 'Value'])
        # Ensure Value column is string type to avoid Arrow serialization issues
        quality_df['Value'] = quality_df['Value'].astype(str)
        st.dataframe(quality_df, width='stretch', hide_index=True)
    
    with tab3:
        st.markdown('<h2 class="section-header">Gene Annotation</h2>', unsafe_allow_html=True)
        
        # Remove duplicates (keep first occurrence) and malformed entries for annotation
        df_clean = df.drop_duplicates(subset=['variant_id'], keep='first').copy()
        # Filter out malformed entries by finding their indices
        duplicates, malformed = validate_variants(df)
        if not malformed.empty:
            malformed_indices = set(malformed['variant_id'].tolist())
            df_clean = df_clean[~df_clean['variant_id'].isin(malformed_indices)]
        
        st.info(f"**{len(df_clean)}** valid variants ready for annotation (duplicates and malformed entries excluded)")
        
        # Initialize session state for annotated data
        if 'annotated_df' not in st.session_state:
            st.session_state.annotated_df = None
        if 'summary' not in st.session_state:
            st.session_state.summary = None
        
        # Annotation button
        col1, col2, col3 = st.columns([1, 1, 1])
        with col2:
            annotate_button = st.button(
                "Annotate Variants with Nearest Genes",
                type="primary",
                width='stretch'
            )
        
        if annotate_button:
            if len(df_clean) == 0:
                st.error("No valid variants to annotate")
            else:
                with st.spinner("Annotating variants... This may take a few minutes."):
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    status_text.info(f"Processing {len(df_clean)} variants. Please wait...")
                    annotated_df = annotate_variants(df_clean, progress_bar)
                    
                    status_text.success("Annotation complete!")
                    progress_bar.empty()
                    
                    # Store in session state
                    st.session_state.annotated_df = annotated_df
                    st.session_state.summary = generate_summary(annotated_df)
                    
                    st.rerun()
        
        # Display annotated data if available
        if st.session_state.annotated_df is not None:
            st.markdown('<h3 style="margin-top: 2rem;">Annotated Variants</h3>', unsafe_allow_html=True)
            
            annotated_pandas = st.session_state.annotated_df.copy()
            
            # Add filters
            st.markdown("**Filter Results**")
            filter_col1, filter_col2, filter_col3 = st.columns(3)
            
            with filter_col1:
                chrom_filter = st.multiselect(
                    "Chromosome",
                    options=sorted(annotated_pandas['chrom'].unique()),
                    default=[]
                )
            
            with filter_col2:
                has_gene = st.selectbox(
                    "Gene Annotation",
                    options=["All", "With Gene", "Without Gene"],
                    index=0
                )
            
            with filter_col3:
                search_variant = st.text_input("Search Variant ID", placeholder="e.g., rs72809442")
            
            # Apply filters
            filtered_df = annotated_pandas.copy()
            
            if chrom_filter:
                filtered_df = filtered_df[filtered_df['chrom'].isin(chrom_filter)]
            
            if has_gene == "With Gene":
                filtered_df = filtered_df[filtered_df['nearest_gene_name'].notna()]
            elif has_gene == "Without Gene":
                filtered_df = filtered_df[filtered_df['nearest_gene_name'].isna()]
            
            if search_variant:
                filtered_df = filtered_df[
                    filtered_df['variant_id'].astype(str).str.contains(search_variant, case=False, na=False)
                ]
            
            st.dataframe(
                filtered_df,
                width='stretch',
                hide_index=True,
                height=500
            )
            
            st.caption(f"Showing {len(filtered_df)} of {len(annotated_pandas)} variants")
    
    with tab4:
        st.markdown('<h2 class="section-header">Summary Statistics</h2>', unsafe_allow_html=True)
        
        if st.session_state.summary is not None:
            summary = st.session_state.summary
            
            # Overall statistics
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Total Variants", summary['total_variants'])
            with col2:
                st.metric("With Annotation", summary['variants_with_gene_annotation'])
            with col3:
                st.metric("Without Annotation", summary['variants_without_gene_annotation'])
            with col4:
                annotation_rate = (summary['variants_with_gene_annotation'] / summary['total_variants'] * 100) if summary['total_variants'] > 0 else 0
                st.metric("Annotation Rate", f"{annotation_rate:.1f}%")
            
            st.divider()
            
            # Chromosome counts with visualization
            st.markdown('<h3>Variants per Chromosome</h3>', unsafe_allow_html=True)
            col1, col2 = st.columns([2, 1])
            
            with col1:
                chrom_counts_dict = summary['chromosome_counts']
                chrom_data = pd.DataFrame({
                    'Chromosome': list(chrom_counts_dict.keys()),
                    'Count': list(chrom_counts_dict.values())
                })
                chrom_data = chrom_data.sort_values('Chromosome')
                
                fig_chrom = px.bar(
                    chrom_data,
                    x='Chromosome',
                    y='Count',
                    title='Variants per Chromosome',
                    labels={'Chromosome': 'Chromosome', 'Count': 'Number of Variants'},
                    color='Count',
                    color_continuous_scale='Viridis'
                )
                fig_chrom.update_layout(
                    showlegend=False,
                    height=400,
                    xaxis_title="Chromosome",
                    yaxis_title="Number of Variants"
                )
                st.plotly_chart(fig_chrom, width='stretch')
            
            with col2:
                st.markdown("**Chromosome Counts**")
                st.dataframe(
                    chrom_data,
                    width='stretch',
                    hide_index=True,
                    height=400
                )
            
            st.divider()
            
            # Top 10 genes with visualization
            st.markdown('<h3>Top 10 Genes with Most Variants</h3>', unsafe_allow_html=True)
            
            if summary['top_10_genes']:
                col1, col2 = st.columns([2, 1])
                
                with col1:
                    gene_counts_dict = summary['top_10_genes']
                    gene_data = pd.DataFrame({
                        'Gene': list(gene_counts_dict.keys()),
                        'Variant Count': list(gene_counts_dict.values())
                    })
                    gene_data = gene_data.sort_values('Variant Count', ascending=True)
                    
                    fig_genes = px.bar(
                        gene_data,
                        x='Variant Count',
                        y='Gene',
                        orientation='h',
                        title='Top 10 Genes by Variant Count',
                        labels={'Variant Count': 'Number of Variants', 'Gene': 'Gene Name'},
                        color='Variant Count',
                        color_continuous_scale='Plasma'
                    )
                    fig_genes.update_layout(
                        showlegend=False,
                        height=400,
                        yaxis={'categoryorder': 'total ascending'}
                    )
                    st.plotly_chart(fig_genes, width='stretch')
                
                with col2:
                    st.markdown("**Gene Rankings**")
                    st.dataframe(
                        gene_data.sort_values('Variant Count', ascending=False),
                        width='stretch',
                        hide_index=True,
                        height=400
                    )
            else:
                st.info("No gene annotations available yet. Please run annotation first.")
            
            st.divider()
            
            # JSON Summary
            st.markdown('<h3>JSON Summary</h3>', unsafe_allow_html=True)
            with st.expander("View/Download JSON Summary", expanded=False):
                st.json(summary)
                
                # Download buttons
                col1, col2 = st.columns(2)
                
                with col1:
                    json_str = json.dumps(summary, indent=2)
                    st.download_button(
                        label="Download Summary JSON",
                        data=json_str,
                        file_name="variant_summary.json",
                        mime="application/json",
                        width='stretch'
                    )
                
                with col2:
                    if st.session_state.annotated_df is not None:
                        csv = st.session_state.annotated_df.to_csv(index=False)
                        st.download_button(
                            label="Download Annotated CSV",
                            data=csv,
                            file_name="annotated_variants.csv",
                            mime="text/csv",
                            width='stretch'
                        )
        else:
            st.info("No summary available yet. Please run annotation first in the 'Annotation' tab.")
else:
    st.error("Failed to load variants from CSV file")
