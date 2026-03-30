import streamlit as st
import pandas as pd
import requests
import plotly.graph_objects as go
import plotly.express as px
import io
import re
import time
import gzip
import numpy as np
import streamlit.components.v1 as components
from streamlit_molstar import st_molstar
import tempfile
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import Polypeptide, is_aa, protein_letters_3to1

st.set_page_config(page_title="Protein Sneakoscope", layout="wide")
st.title("🧬 The Protein Sneakoscope")

# --- 1. INITIALIZE SESSION STATE ---
# This acts as our "memory" so data doesn't disappear
if 'data_loaded' not in st.session_state:
    st.session_state.data_loaded = False
    st.session_state.var_df = None
    st.session_state.am_df = None
    st.session_state.feat_df = None
    st.session_state.seq_len = 0
    st.session_state.sequence = ""
    st.session_state.uniprot_id = ""

# --- 2. CACHED DATA FETCHING ---
@st.cache_data
def fetch_all_data(uid):
    """Fetches all data in one go and caches it."""
    # 1. Features
    feat_url = f"https://www.ebi.ac.uk/proteins/api/features/{uid}"
    feat_res = requests.get(feat_url, headers={"Accept": "application/json"})
    feat_data, s_len, sequence = None, 0, ""
    if feat_res.status_code == 200:
        d = feat_res.json()
        s_len = int(d.get('sequenceLength', 0))
        sequence = d.get('sequence', '')
        feat_data = pd.DataFrame([{'type': f.get('type'), 'description': f.get('description', 'N/A'),
                                   'start': int(f.get('begin')), 'end': int(f.get('end'))} 
                                  for f in d.get('features', [])])
    
    # 2. Variants
    var_url = f"https://www.ebi.ac.uk/proteins/api/variation/{uid}"
    var_res = requests.get(var_url, headers={"Accept": "application/json"})
    var_data = None
    if var_res.status_code == 200:
        v_feats = var_res.json().get('features', [])
        recs = []
        for f in v_feats:
            w, m, p = f.get('wildType'), f.get('alternativeSequence'), f.get('begin')
            if w and m and p and len(w) == 1 and len(m) == 1:
                recs.append({'protein_variant': f"{w}{p}{m}", 'position': int(p), 'source': f.get('sourceType')})
        if recs: var_data = pd.DataFrame(recs)

    # 3. AlphaMissense
    am_data = None
    for url in [f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-aa_substitutions.csv",
                f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-hg38.csv"]:
        r = requests.get(url)
        if r.status_code == 200:
            am_data = pd.read_csv(io.StringIO(r.text))
            break
            
    return var_data, am_data, feat_data, s_len, sequence

# --- 3.5. RAMACHANDRAN TOOLS ---
def get_color_for_plddt(plddt):
    if plddt > 90: return "#0053D6"
    elif plddt > 70: return "#65CBF3"
    elif plddt > 50: return "#FFDB13"
    else: return "#FF7D45"

def get_ramachandran_df(pdb_text, text_id):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(text_id, io.StringIO(pdb_text))
    except Exception as e:
        return pd.DataFrame()
        
    data = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    for atom in residue:
                        if atom.get_name() == 'CA':
                            data.append({
                                'Residue': residue.id[1],
                                'pLDDT': atom.get_bfactor(),
                                'AA': protein_letters_3to1.get(residue.get_resname().upper(), 'X')
                            })
                            break
                            
    if not data: return pd.DataFrame()
    df = pd.DataFrame(data).sort_values(by='Residue').reset_index(drop=True)
    
    phi_psi_list = []
    for model in structure:
        for chain in model:
            residues = [res for res in chain if is_aa(res, standard=True)]
            if len(residues) < 2: continue
            poly_chain = Polypeptide(residues)
            phi_psi_tuples = poly_chain.get_phi_psi_list()
            for i, res in enumerate(poly_chain):
                res_id = res.get_id()[1]
                phi, psi = phi_psi_tuples[i]
                phi_psi_list.append({
                    'Residue': res_id,
                    'Phi': np.degrees(phi) if phi is not None else None,
                    'Psi': np.degrees(psi) if psi is not None else None,
                })
                
    if not phi_psi_list:
        df[['Phi', 'Psi']] = None
        return df
        
    phi_psi_df = pd.DataFrame(phi_psi_list)
    return pd.merge(df, phi_psi_df, on='Residue', how='left')

def generate_ramachandran_plot(df, title):
    plot_df = df.dropna(subset=['Phi', 'Psi']).copy()
    if plot_df.empty: return None
    plot_df['Color'] = plot_df['pLDDT'].apply(get_color_for_plddt)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=plot_df['Phi'], y=plot_df['Psi'], mode='markers',
        marker=dict(color=plot_df['Color'], size=6, line=dict(width=0.5, color='DarkSlateGrey'), showscale=False),
        text=plot_df.apply(lambda r: f"Residue: {r['AA']}{r['Residue']}<br>pLDDT: {r['pLDDT']:.2f}<br>Phi: {r['Phi']:.2f}<br>Psi: {r['Psi']:.2f}", axis=1),
        hoverinfo='text', name='Residues'
    ))
    
    shapes = [
        dict(type="rect", xref="x", yref="y", x0=-180, y0=100, x1=-40, y1=180, fillcolor="rgba(173, 216, 230, 0.2)", layer="below", line_width=0),
        dict(type="rect", xref="x", yref="y", x0=-160, y0=-70, x1=-30, y1=50, fillcolor="rgba(144, 238, 144, 0.2)", layer="below", line_width=0),
        dict(type="rect", xref="x", yref="y", x0=30, y0=0, x1=100, y1=100, fillcolor="rgba(255, 182, 193, 0.2)", layer="below", line_width=0)
    ]
    
    fig.update_layout(
        title=title, xaxis_title="Phi (Φ) degrees", yaxis_title="Psi (Ψ) degrees",
        xaxis=dict(range=[-180, 180], tickvals=[-180, -120, -60, 0, 60, 120, 180], zeroline=True, zerolinecolor='black', zerolinewidth=1),
        yaxis=dict(range=[-180, 180], tickvals=[-180, -120, -60, 0, 60, 120, 180], zeroline=True, zerolinecolor='black', zerolinewidth=1),
        height=450, showlegend=False, shapes=shapes, margin=dict(l=40, r=40, t=40, b=40), template='plotly_white'
    )
    return fig

# --- 3. SIDEBAR CONTROLS ---
st.sidebar.header("📜 Controls")
uniprot_id = st.sidebar.text_input("Enter UniProt Accession", value="P63000").strip().upper()

st.sidebar.markdown("---")
st.sidebar.subheader("🛠️ SWISS-MODEL Integration")
st.session_state.swiss_token = st.sidebar.text_input("API Token (Optional)", type="password", help="Get a free token from your SWISS-MODEL account profile.")

if st.sidebar.button("Cast Revelio!"):
    with st.spinner("Consulting the Great Library..."):
        v, a, f, s, seq = fetch_all_data(uniprot_id)
        st.session_state.var_df = v
        st.session_state.am_df = a
        st.session_state.feat_df = f
        st.session_state.seq_len = s
        st.session_state.sequence = seq
        st.session_state.uniprot_id = uniprot_id
        st.session_state.data_loaded = True

# --- 4. MAIN DISPLAY LOGIC ---
# This runs every time the page refreshes, as long as data_loaded is True
if st.session_state.data_loaded:
    var_df = st.session_state.var_df
    feat_df = st.session_state.feat_df
    am_df = st.session_state.am_df
    seq_len = st.session_state.seq_len

    # A. Mutation Frequency Plot
    st.subheader("📊 Mutation Frequency per Site")
    freq_df = var_df['position'].value_counts().reset_index()
    freq_df.columns = ['position', 'count']
    fig_freq = px.bar(freq_df.sort_values('position'), x='position', y='count', color_discrete_sequence=['#636EFA'])
    fig_freq.update_layout(height=250, margin=dict(t=10, b=10))
    st.plotly_chart(fig_freq, use_container_width=True)

# B. Interactive Feature Map (Fixed Y-Axis and Range)
    if feat_df is not None:
        st.subheader("🗺️ Feature Map: Site Architecture")
        all_types = sorted(feat_df['type'].unique().tolist())
        selected_types = st.sidebar.multiselect("Select Features to Reveal:", all_types, default=all_types[:3])
        
        if selected_types:
            fig_map = go.Figure()
            
            # 1. Sequence Backbone at y=0
            fig_map.add_trace(go.Scatter(
                x=[1, seq_len], y=[0, 0], 
                mode="lines", 
                line=dict(color="black", width=2), 
                hoverinfo="skip", showlegend=False
            ))
            
            color_palette = px.colors.qualitative.Prism
            for i, f_type in enumerate(selected_types):
                type_subset = feat_df[feat_df['type'] == f_type]
                color = color_palette[i % len(color_palette)]
                
                # Stagger traces: first feature at y=1, second at y=2, etc.
                y_val = i + 1 
                
                for _, row in type_subset.iterrows():
                    fig_map.add_trace(go.Scatter(
                        x=[row['start'], row['end']], 
                        y=[y_val, y_val],
                        mode="lines", 
                        line=dict(color=color, width=15),
                        name=f_type, 
                        text=f"<b>{f_type}</b><br>{row['description']}<br>Pos: {row['start']}-{row['end']}",
                        hoverinfo="text", 
                        legendgroup=f_type, 
                        showlegend=True if _ == type_subset.index[0] else False
                    ))
            
            # 2. DYNAMIC RANGE FIX
            # We set the range from -0.5 to (Number of Features + 1)
            # This ensures every trace has its own visible lane.
            fig_map.update_layout(
                height=200 + (len(selected_types) * 40), # Plot gets taller as you add features
                xaxis_title="Residue Position",
                yaxis=dict(
                    showticklabels=False, 
                    range=[-0.5, len(selected_types) + 1], # Ensures all y_vals are visible
                    fixedrange=True # Prevents accidental vertical zooming
                ),
                margin=dict(l=20, r=20, t=10, b=10),
                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                hovermode="closest"
            )
            st.plotly_chart(fig_map, use_container_width=True)
        else:
            st.warning("Please select at least one feature in the sidebar to reveal the map.")

    # C. Pathogenicity Plot
    if am_df is not None:
        st.subheader("🎯 Pathogenicity Analysis")
        merged_df = pd.merge(var_df, am_df, on='protein_variant', how='inner')
        fig_patho = px.scatter(merged_df, x='position', y='am_pathogenicity', color='am_class',
                               color_discrete_map={'likely_benign': 'green', 'ambiguous': 'grey', 'likely_pathogenic': 'red'})
        fig_patho.add_hline(y=0.34, line_dash="dash", line_color="green", opacity=0.3)
        fig_patho.add_hline(y=0.564, line_dash="dash", line_color="red", opacity=0.3)
        st.plotly_chart(fig_patho, use_container_width=True)

    # D. Mutation Selection and Alignment
    st.subheader("✂️ Mutation Constructor & Alignment")
    ref_seq = st.session_state.sequence
    if ref_seq and var_df is not None and not var_df.empty:
        available_variants = sorted(var_df['protein_variant'].dropna().unique())
        selected_variants = st.multiselect("Select Mutations to Apply", available_variants)
        
        mutated_seq_list = list(ref_seq)
        mutations_applied = []
        for var in selected_variants:
            match = re.match(r"([a-zA-Z]+)(\d+)([a-zA-Z*]+)", var)
            if match:
                wt, pos, mut = match.groups()
                pos_idx = int(pos) - 1
                if 0 <= pos_idx < len(mutated_seq_list) and mutated_seq_list[pos_idx] == wt:
                    mutated_seq_list[pos_idx] = mut
                    mutations_applied.append(var)
                else:
                    st.warning(f"Mutation {var} does not match reference {wt} at position {pos}")
        
        mutated_seq = "".join(mutated_seq_list)
        
        if mutations_applied:
            st.markdown("### 🔮 AlphaMissense Consequence Predictions")
            am_df = st.session_state.get('am_df')
            if am_df is not None and isinstance(am_df, pd.DataFrame) and not am_df.empty:
                for var in mutations_applied:
                    match_row = am_df[am_df['protein_variant'] == var]
                    if not match_row.empty:
                        am_class = str(match_row.iloc[0]['am_class'])
                        am_score = float(match_row.iloc[0]['am_pathogenicity'])
                        color = "#008000" if am_class == 'likely_benign' else "#FF0000" if am_class == 'likely_pathogenic' else "#808080"
                        st.markdown(f"- **{var}**: <span style='color:{color}; font-weight:bold;'>{am_class.replace('_', ' ').title()}</span> *(Score: {am_score:.3f})*", unsafe_allow_html=True)
                    else:
                        st.markdown(f"- **{var}**: *Data not available in AlphaMissense*")
            else:
                st.info("AlphaMissense data is unavailable for these variants.")
                
            st.markdown(f"**Mutated FASTA ({len(mutations_applied)} mutations applied)**")
            fasta_header = f">{st.session_state.uniprot_id} | Mutated | {', '.join(mutations_applied)}"
            fasta_output = f"{fasta_header}\n{mutated_seq}"
            st.code(fasta_output, language='fasta')
            
            st.markdown("**Visualized Alignment**")
            html = "<div style='font-family: monospace; font-size: 14px; white-space: pre-wrap; background-color: #f8f9fa; padding: 15px; border-radius: 5px; line-height: 1.5; color: black'>"
            chunk_size = 60
            for i in range(0, len(ref_seq), chunk_size):
                ref_chunk = ref_seq[i:i+chunk_size]
                mut_chunk = mutated_seq[i:i+chunk_size]
                
                html += f"<b>Ref {i+1:04d}:</b> "
                for r_aa in ref_chunk:
                    html += r_aa
                html += "<br>"
                
                html += f"<b>Mut {i+1:04d}:</b> "
                for j, (r_aa, m_aa) in enumerate(zip(ref_chunk, mut_chunk)):
                    if r_aa != m_aa:
                        html += f"<span style='background-color: #ffcccc; color: #cc0000; font-weight: bold;'>{m_aa}</span>"
                    else:
                        html += m_aa
                html += "<br><br>"
                
            html += "</div>"
            st.markdown(html, unsafe_allow_html=True)
            
            # --- SWISS-MODEL INTEGRATION ---
            st.markdown("---")
            st.subheader("🧬 3D Structure Prediction")
            if st.button("Generate 3D Model with SWISS-MODEL"):
                if not st.session_state.swiss_token:
                    st.error("Please enter your SWISS-MODEL API Token in the sidebar first.")
                else:
                    with st.spinner("Downloading AlphaFold Template..."):
                        af_pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{st.session_state.uniprot_id}-F1-model_v6.pdb"
                        af_req = requests.get(af_pdb_url)
                        if af_req.status_code != 200:
                            st.error(f"Could not retrieve AlphaFold structure for {st.session_state.uniprot_id}")
                            st.stop()
                        st.session_state.af_pdb_text = af_req.text
                    
                    with st.spinner("Submitting job to SWISS-MODEL... this will take a few minutes!"):
                        sm_url = "https://swissmodel.expasy.org/user_template/"
                        headers = {
                            "Authorization": f"Token {st.session_state.swiss_token}",
                            "Content-Type": "application/json"
                        }
                        payload = {
                            "target_sequences": [mutated_seq],
                            "project_title": f"{st.session_state.uniprot_id} Mutant",
                            "template_coordinates": st.session_state.af_pdb_text
                        }
                        
                        try:
                            # 1. Submit
                            res = requests.post(sm_url, headers=headers, json=payload)
                            if res.status_code not in (200, 201, 202):
                                st.error(f"Failed to submit to SWISS-MODEL: {res.text}")
                                st.stop()
                            
                            proj_data = res.json()
                            status_url = proj_data.get('url')
                            if not status_url and 'project_id' in proj_data:
                                status_url = f"https://swissmodel.expasy.org/project/{proj_data['project_id']}/models/summary/"
                                
                            if not status_url:
                                st.error("No project URL returned from SWISS-MODEL.")
                                st.stop()
                                
                            # 2. Poll Status
                            completed = False
                            for step in range(60): # wait up to 10 mins (60 * 10s)
                                time.sleep(10)
                                stat_res = requests.get(status_url, headers=headers)
                                if stat_res.status_code != 200:
                                    continue
                                stat_data = stat_res.json()
                                status = stat_data.get('status', '')
                                if status == 'COMPLETED':
                                    completed = True
                                    results = stat_data
                                    break
                                elif status in ['FAILED', 'CANCELLED', 'ERROR']:
                                    st.error(f"SWISS-MODEL Job Failed: {stat_data.get('failure_message', status)}")
                                    st.stop()
                                    
                            if not completed:
                                st.error("Timeout waiting for SWISS-MODEL. Please check your dashboard on their website.")
                                st.stop()
                                
                            # 3. Retrieve results
                            models = results.get('models', [])
                            if not models:
                                st.error("Job completed but no models were built.")
                                st.stop()
                                
                            model_coord_url = models[0].get('coordinates_url')
                            if model_coord_url:
                                pdb_res = requests.get(model_coord_url, headers=headers)
                                if pdb_res.content.startswith(b'\x1f\x8b'):
                                    st.session_state.mut_pdb_text = gzip.decompress(pdb_res.content).decode("utf-8")
                                else:
                                    st.session_state.mut_pdb_text = pdb_res.content.decode("utf-8")
                                
                                st.session_state.sm_stats = models[0]
                                st.success("3D Model successfully generated!")
                                
                        except Exception as e:
                            st.error(f"Error communicating with SWISS-MODEL: {str(e)}")
                            
            # Display 3D comparison if available
            if 'mut_pdb_text' in st.session_state and 'af_pdb_text' in st.session_state:
                st.markdown("### 🔍 Interactive 3D Comparison")
                st.markdown("**(Left: AlphaFold Wild-Type | Right: Mutated SWISS-MODEL)**")
                st.info("Colored by pLDDT (B-factor). Red: low confidence, Blue: high confidence. Rotate with your mouse to sync both views.")
                
                # Download button
                st.download_button("Download Mutated PDB", data=st.session_state.mut_pdb_text, file_name=f"{st.session_state.uniprot_id}_mutant.pdb", mime="text/plain")
                
                st.markdown("---")
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**AlphaFold Wild-Type**")
                    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp1:
                        tmp1.write(st.session_state.af_pdb_text.encode('utf-8'))
                        tmp1_path = tmp1.name
                    st_molstar(tmp1_path, height='500px', key='wt_molstar')
                    
                with col2:
                    st.markdown("**Mutated SWISS-MODEL**")
                    if 'sm_stats' in st.session_state:
                         qmean = st.session_state.sm_stats.get('qmean', 'N/A')
                         gmqe = st.session_state.sm_stats.get('gmqe', 'N/A')
                         if qmean != 'N/A': st.caption(f"**SWISS-MODEL Est:** QMEAN: {qmean:.3f} | GMQE: {gmqe:.3f}")
                         
                    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp2:
                        tmp2.write(st.session_state.mut_pdb_text.encode('utf-8'))
                        tmp2_path = tmp2.name
                    st_molstar(tmp2_path, height='500px', key='mut_molstar')

                # --- RAMACHANDRAN PLOTS ---
                st.markdown("### 📈 Structure Quality Assessment (Ramachandran)")
                rc1, rc2 = st.columns(2)
                with rc1:
                    df_wt = get_ramachandran_df(st.session_state.af_pdb_text, "WT")
                    fig_wt = generate_ramachandran_plot(df_wt, "AlphaFold Wild-Type")
                    if fig_wt: st.plotly_chart(fig_wt, use_container_width=True)
                with rc2:
                    df_mut = get_ramachandran_df(st.session_state.mut_pdb_text, "Mutant")
                    fig_mut = generate_ramachandran_plot(df_mut, "Mutated SWISS-MODEL")
                    if fig_mut: st.plotly_chart(fig_mut, use_container_width=True)

        else:
            st.info("Select one or more mutations from the list above to generate the custom sequence and alignment.")
    else:
        st.warning("Sequence or variant data is unavailable for this protein.")

else:
    st.info("Enter a UniProt ID and click 'Cast Revelio!' to start.")