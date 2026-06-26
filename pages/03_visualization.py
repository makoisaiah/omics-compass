# OmicsCompass - 可視化ページ
# Plotly を使ったインタラクティブな可視化:
# https://plotly.com/python/
# Volcano plot の作り方:
# https://plotly.com/python/volcano-plot/
# gseapy を使ったパスウェイ解析:
# https://gseapy.readthedocs.io/en/latest/run.html

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import gseapy as gp
import networkx as nx
import requests

st.title("📊 Visualization")
st.subheader("パスウェイ解析と可視化")

# セッションにデータがあるか確認
if "df_results" not in st.session_state:
    st.warning("まず Analysis ページで差分解析を実行してください")
    st.stop()

df_results = st.session_state["df_results"]
df_sig = st.session_state["df_sig"]

st.markdown(f"**有意な変化のあったプローブ数**: {len(df_sig)}")

# タブで可視化を切り替え
tab1, tab2, tab3, tab4 = st.tabs(["🌋 Volcano Plot", "🔥 Heatmap", "🗺️ Pathway Analysis", "🕸️ Network"])

# --- タブ1: Volcano Plot ---
with tab1:
    st.markdown("### Volcano Plot")
    st.markdown("横軸: Log2 Fold Change、縦軸: -log10(p値)")

    df_plot = df_results.copy()
    df_plot["-log10(p_adj)"] = -np.log10(df_plot["p_adj"].clip(lower=1e-300))

    # 色分け
    def classify(row):
        if row["p_adj"] < 0.05 and row["Log2FoldChange"] > 1:
            return "上昇 (Up)"
        elif row["p_adj"] < 0.05 and row["Log2FoldChange"] < -1:
            return "低下 (Down)"
        else:
            return "変化なし"

    df_plot["分類"] = df_plot.apply(classify, axis=1)

    color_map = {
        "上昇 (Up)": "#e74c3c",
        "低下 (Down)": "#3498db",
        "変化なし": "#95a5a6"
    }

    fig = px.scatter(
        df_plot,
        x="Log2FoldChange",
        y="-log10(p_adj)",
        color="分類",
        color_discrete_map=color_map,
        hover_data=["Gene", "Probe", "p_adj"],
        title="Volcano Plot",
        labels={
            "Log2FoldChange": "Log2 Fold Change",
            "-log10(p_adj)": "-log10(Adjusted p-value)"
        }
    )

    # 閾値ラインを追加
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray",
                  annotation_text="p_adj = 0.05")
    fig.add_vline(x=1, line_dash="dash", line_color="gray")
    fig.add_vline(x=-1, line_dash="dash", line_color="gray")

    st.plotly_chart(fig, use_container_width=True)

# --- タブ2: Heatmap ---
with tab2:
    st.markdown("### Heatmap")
    st.markdown("有意に変化した上位プローブの発現パターン")

    if "gse" not in st.session_state:
        st.warning("Data Fetch ページでデータを再取得してください")
    else:
        gse = st.session_state["gse"]
        top_probes = df_sig.head(30)["Probe"].tolist()
        probe_to_gene = df_sig.set_index("Probe")["Gene"].to_dict()

        # 発現データを再取得
        expr_data = {}
        for name, gsm in gse.gsms.items():
            if gsm.table is not None and not gsm.table.empty:
                table = gsm.table.set_index("ID_REF")["VALUE"]
                expr_data[name] = table

        df_expr = pd.DataFrame(expr_data)
        df_expr = df_expr.apply(pd.to_numeric, errors="coerce").dropna()

        # 上位プローブだけ抽出
        df_heat = df_expr.loc[
            df_expr.index.isin(top_probes)
        ].copy()

        if not df_heat.empty:
            # log2 変換
            df_heat = np.log2(df_heat.clip(lower=1e-10))
            df_heat.index = df_heat.index.map(lambda x: probe_to_gene.get(x, x))

            fig2 = px.imshow(
                df_heat,
                aspect="auto",
                color_continuous_scale="RdBu_r",
                title="発現ヒートマップ（上位30プローブ）",
                labels={"color": "log2 発現量"}
            )
            st.plotly_chart(fig2, use_container_width=True)

# --- タブ3: Pathway Analysis ---
with tab3:
    st.markdown("### パスウェイ解析")
    st.markdown("有意に変化したプローブを使って関連パスウェイを探索します")

    # プローブ ID を遺伝子シンボルに変換する必要があるが
    # まずプローブリストで試す
    gene_list = df_sig["Gene"].tolist()

    db_options = [
        "KEGG_2021_Human",
        "Reactome_2022",
        "GO_Biological_Process_2021",
        "WikiPathway_2021_Human"
    ]

    selected_db = st.selectbox("データベースを選択", db_options)

    if st.button("パスウェイ解析を実行", type="primary"):
        with st.spinner("解析中..."):
            try:
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=selected_db,
                    organism="mouse",  # GSE21393 はマウスのデータ
                    outdir=None
                )

                df_enr = enr.results.sort_values("Adjusted P-value")
                df_enr_sig = df_enr[df_enr["Adjusted P-value"] < 0.05]

                if df_enr_sig.empty:
                    st.warning("有意なパスウェイが見つかりませんでした。プローブIDから遺伝子シンボルへの変換が必要かもしれません。")
                else:
                    st.success(f"{len(df_enr_sig)} 個の有意なパスウェイが見つかりました")

                    # 上位パスウェイをバーチャートで表示
                    df_top = df_enr_sig.head(15).copy()
                    df_top["-log10(p_adj)"] = -np.log10(
                        df_top["Adjusted P-value"].clip(lower=1e-300)
                    )

                    fig3 = px.bar(
                        df_top.sort_values("-log10(p_adj)"),
                        x="-log10(p_adj)",
                        y="Term",
                        orientation="h",
                        title=f"上位パスウェイ ({selected_db})",
                        color="-log10(p_adj)",
                        color_continuous_scale="Reds"
                    )
                    st.plotly_chart(fig3, use_container_width=True)
                    st.dataframe(df_enr_sig[["Term", "Adjusted P-value", "Genes"]].head(20))

            except Exception as e:
                st.error(f"エラーが発生しました: {e}")
                st.exception(e)

# --- タブ4: ネットワーク解析 ---
with tab4:
    st.markdown("### 遺伝子ネットワーク解析")
    st.markdown("""
    STRING DB を使って、有意に変化した遺伝子間の相互作用ネットワークを表示します。
    - **線が太い**: 相互作用の信頼度が高い
    - **赤いノード**: 発現上昇
    - **青いノード**: 発現低下
    """)

    # STRING DB API について:
    # https://string-db.org/help/api/

    # 上位遺伝子を取得
    n_genes = st.slider("解析する遺伝子数", min_value=10, max_value=50, value=20)
    
    top_genes = df_sig.head(n_genes)
    gene_symbols = top_genes["Gene"].tolist()
    
    # プローブIDが混入している場合は除外
    gene_symbols = [g for g in gene_symbols if not g.endswith("_at")]

    if st.button("ネットワーク解析を実行", type="primary"):
        with st.spinner("STRING DB に問い合わせ中..."):
            try:
                # STRING DB API でネットワークデータを取得
                string_api_url = "https://string-db.org/api/json/network"
                params = {
                    "identifiers": "%0d".join(gene_symbols),
                    "species": 10090,  # マウスの taxonomy ID
                    "caller_identity": "omics_compass"
                }
                response = requests.post(string_api_url, data=params)
                data = response.json()

                if not data:
                    st.warning("ネットワークデータが見つかりませんでした。遺伝子シンボルの変換を確認してください。")
                    st.stop()

                # NetworkX でグラフを構築
                G = nx.Graph()

                # ノードを追加（発現変化の情報付き）
                gene_fc = df_sig.set_index("Gene")["Log2FoldChange"].to_dict()
                for gene in gene_symbols:
                    G.add_node(gene, log2fc=gene_fc.get(gene, 0))

                # エッジを追加
                for interaction in data:
                    gene_a = interaction["preferredName_A"]
                    gene_b = interaction["preferredName_B"]
                    score = interaction["score"]
                    if score > 0.4:  # 信頼度 0.4 以上のみ
                        G.add_edge(gene_a, gene_b, weight=score)

                st.success(f"ネットワーク構築完了: {G.number_of_nodes()} 遺伝子, {G.number_of_edges()} 相互作用")

                # Plotly でネットワークを可視化
                # spring layout でノードの位置を決定
                pos = nx.spring_layout(G, seed=42)

                # エッジの描画
                edge_traces = []
                for edge in G.edges(data=True):
                    x0, y0 = pos[edge[0]]
                    x1, y1 = pos[edge[1]]
                    weight = edge[2].get("weight", 0.5)
                    edge_traces.append(
                        go.Scatter(
                            x=[x0, x1, None],
                            y=[y0, y1, None],
                            mode="lines",
                            line=dict(width=weight * 3, color="#888"),
                            hoverinfo="none",
                            showlegend=False
                        )
                    )

                # ノードの描画
                node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
                for node in G.nodes():
                    x, y = pos[node]
                    node_x.append(x)
                    node_y.append(y)
                    fc = G.nodes[node].get("log2fc", 0)
                    node_text.append(f"{node}<br>Log2FC: {fc:.2f}")
                    node_color.append(fc)
                    # 次数（接続数）でサイズを変える
                    node_size.append(10 + G.degree(node) * 5)

                node_trace = go.Scatter(
                    x=node_x,
                    y=node_y,
                    mode="markers+text",
                    text=list(G.nodes()),
                    textposition="top center",
                    hovertext=node_text,
                    hoverinfo="text",
                    marker=dict(
                        size=node_size,
                        color=node_color,
                        colorscale="RdBu_r",
                        colorbar=dict(title="Log2FC"),
                        line=dict(width=1, color="white")
                    )
                )

                fig4 = go.Figure(
                    data=edge_traces + [node_trace],
                    layout=go.Layout(
                        title="遺伝子相互作用ネットワーク（STRING DB）",
                        showlegend=False,
                        hovermode="closest",
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        height=600
                    )
                )

                st.plotly_chart(fig4, use_container_width=True)

                # ハブ遺伝子（接続数が多い遺伝子）を表示
                st.markdown("### ハブ遺伝子（相互作用が多い遺伝子）")
                degree_df = pd.DataFrame(
                    [(node, G.degree(node), gene_fc.get(node, 0))
                     for node in G.nodes()],
                    columns=["Gene", "接続数", "Log2FoldChange"]
                ).sort_values("接続数", ascending=False)

                st.dataframe(degree_df.head(10), use_container_width=True)

            except Exception as e:
                st.error(f"エラーが発生しました: {e}")
                st.exception(e)