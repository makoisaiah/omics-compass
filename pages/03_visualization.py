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
tab1, tab2, tab3 = st.tabs(["🌋 Volcano Plot", "🔥 Heatmap", "🗺️ Pathway Analysis"])

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
        hover_data=["Probe", "p_adj"],
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
    gene_list = df_sig["Probe"].tolist()

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