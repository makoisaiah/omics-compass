# OmicsCompass - 差分解析ページ
# マイクロアレイデータの差分解析について:
# https://geoparse.readthedocs.io/en/latest/
# scipy を使った t 検定:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
# 多重検定補正 (Benjamini-Hochberg):
# https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
# mygene を使ったプローブ変換:
# https://docs.mygene.info/projects/mygene-py/en/latest/

import streamlit as st
import pandas as pd
import numpy as np
import json
import re
from scipy import stats
from statsmodels.stats.multitest import multipletests
import mygene
import requests
from components.llm import ask_llm, is_ollama_available

st.title("🔬 Analysis")
st.subheader("差分発現解析")

if "gse" not in st.session_state:
    st.warning("まず Data Fetch ページでデータを取得してください")
    st.stop()

gse = st.session_state["gse"]
geo_id = st.session_state.get("geo_id", "")

st.markdown(f"**データセット**: {geo_id}")

sample_names = list(gse.gsms.keys())

# df_expr がセッションにある場合、実際に存在するサンプルのみに絞り込む
if "df_expr" in st.session_state:
    available = st.session_state["df_expr"].columns.tolist()
    sample_names = [n for n in sample_names if n in available]
    if not sample_names:
        # カラム名が GSM ID でない場合（例: DMSO_HDM_rep1）は全サンプルを使う
        sample_names = list(gse.gsms.keys())

sample_titles = {
    name: gse.gsms[name].metadata.get("title", ["不明"])[0]
    for name in sample_names
}

# LLM によるグループ自動推定
st.markdown("### グループ自動推定")
backend = "ローカル Mistral" if is_ollama_available() else "Groq API"
st.info(f"使用する LLM: {backend}")

if st.button("LLM でグループを自動判定"):
    sample_list = "\n".join([
        f"{name}: {sample_titles[name]}"
        for name in sample_names
    ])
    prompt = f"""You are a bioinformatics expert. Classify each RNA-seq sample as either control or treatment.

Rules:
- vehicle, DMSO, PBS, wild type, WT, sham, untreated, control → "コントロール"
- drug, inhibitor, knockdown, knockout, KO, overexpression, treated, mutant → "処理群"
- If unclear, use the biological context to decide

Sample list:
{sample_list}

Reply ONLY with valid JSON. No explanation. No markdown. Example:
{{"GSM001": "コントロール", "GSM002": "処理群"}}

Your answer:"""

    with st.spinner("LLM が判定中..."):
        try:
            groq_key = None
            try:
                groq_key = st.secrets["GROQ_API_KEY"]
            except Exception:
                pass
            result, used_backend = ask_llm(prompt, groq_key)
            st.success(f"{used_backend} で判定しました")

            json_match = re.search(r'\{.*\}', result, re.DOTALL)
            if json_match:
                suggestions = json.loads(json_match.group())
                st.session_state["group_suggestions"] = suggestions
                for sample_name, group_value in suggestions.items():
                    st.session_state[f"group_{sample_name}"] = group_value
                st.markdown("**判定結果（確認して必要なら修正してください）**")
                st.json(suggestions)
            else:
                st.warning("JSON の解析に失敗しました。手動で設定してください。")
        except Exception as e:
            st.error(f"エラー: {e}")

st.divider()

# グループ割り当て UI
st.markdown("### グループ分け")
group_assignments = {}
suggestions = st.session_state.get("group_suggestions", {})

cols = st.columns(2)
with cols[0]:
    st.markdown("**コントロール群**")
with cols[1]:
    st.markdown("**処理群**")

for name in sample_names:
    title = sample_titles[name]
    col1, col2 = st.columns(2)
    with col1:
        st.text(f"{name}: {title[:40]}")
    with col2:
        suggested = suggestions.get(name, "コントロール")
        if f"group_{name}" not in st.session_state:
            st.session_state[f"group_{name}"] = suggested
        group = st.selectbox(
            label=name,
            options=["コントロール", "処理群"],
            key=f"group_{name}",
            label_visibility="collapsed"
        )
        group_assignments[name] = group

if suggestions:
    n_control = sum(1 for v in suggestions.values() if v == "コントロール")
    n_treatment = sum(1 for v in suggestions.values() if v == "処理群")
    st.caption(f"LLM 判定: コントロール {n_control} サンプル、処理群 {n_treatment} サンプル　※必要に応じて手動で修正できます")

# 解析実行ボタン
if st.button("差分解析を実行", type="primary"):
    control_samples = [n for n, g in group_assignments.items() if g == "コントロール"]
    treatment_samples = [n for n, g in group_assignments.items() if g == "処理群"]

    if len(control_samples) < 2 or len(treatment_samples) < 2:
        st.error("各群に最低2サンプル必要です")
        st.stop()

    with st.spinner("解析中..."):
        try:
            # 発現データを取得
            if "df_expr" in st.session_state:
                df_expr = st.session_state["df_expr"].copy()
                # カラム名が GSM ID でない場合、順番で対応付け
                if not any(s in df_expr.columns for s in control_samples):
                    col_map = dict(zip(df_expr.columns, sample_names))
                    df_expr = df_expr.rename(columns=col_map)
                    st.info(f"カラム対応付け: {col_map}")
            else:
                # マイクロアレイの場合は GSM テーブルから取得
                expr_data = {}
                for name in sample_names:
                    gsm = gse.gsms[name]
                    if gsm.table is not None and not gsm.table.empty:
                        table = gsm.table.set_index("ID_REF")["VALUE"]
                        expr_data[name] = table

                if not expr_data:
                    st.error("発現データが見つかりません。Data Fetch ページで「発現データを取得」ボタンを押してください。")
                    st.stop()

                df_expr = pd.DataFrame(expr_data).dropna()
                df_expr = df_expr.apply(pd.to_numeric, errors="coerce").dropna()

            st.success(f"発現データ取得完了: {len(df_expr)} プローブ")

            # 利用可能なサンプルに絞り込む
            available_samples = df_expr.columns.tolist()
            control_samples = [s for s in control_samples if s in available_samples]
            treatment_samples = [s for s in treatment_samples if s in available_samples]

            if not control_samples or not treatment_samples:
                st.error(f"有効なサンプルが不足しています。利用可能: {available_samples}")
                st.stop()

            st.info(f"解析に使用するサンプル → コントロール: {control_samples} / 処理群: {treatment_samples}")

            # t 検定で差分解析
            control_data = df_expr[control_samples].values
            treatment_data = df_expr[treatment_samples].values

            t_stats, p_values = stats.ttest_ind(
                treatment_data, control_data, axis=1
            )

            # Fold Change を計算
            mean_control = np.mean(control_data, axis=1)
            mean_treatment = np.mean(treatment_data, axis=1)
            mean_control_safe = np.where(mean_control > 0, mean_control, 1e-10)
            mean_treatment_safe = np.where(mean_treatment > 0, mean_treatment, 1e-10)
            log2fc = np.log2(mean_treatment_safe) - np.log2(mean_control_safe)

            # 多重検定補正
            _, p_adj, _, _ = multipletests(p_values, method="fdr_bh")

            df_results = pd.DataFrame({
                "Probe": df_expr.index,
                "Log2FoldChange": log2fc,
                "p_value": p_values,
                "p_adj": p_adj,
            }).sort_values("p_adj")

            # サンプル数が少ない場合は閾値を緩める
            n_min = min(len(control_samples), len(treatment_samples))
            if n_min <= 2:
                p_threshold = 0.1
                fc_threshold = 0.5
                st.info(f"サンプル数が少ないため閾値を緩めています（p_adj < {p_threshold}, |Log2FC| > {fc_threshold}）")
            else:
                p_threshold = 0.05
                fc_threshold = 1.0

            df_sig = df_results[
                (df_results["p_adj"] < p_threshold) &
                (abs(df_results["Log2FoldChange"]) > fc_threshold)
            ].copy()

            # プローブ ID → 遺伝子シンボル変換
            with st.spinner("プローブ ID を遺伝子シンボルに変換中..."):
                mg = mygene.MyGeneInfo()
                probe_ids = df_results["Probe"].tolist()

                result = mg.querymany(
                    probe_ids,
                    scopes="reporter,symbol",
                    fields="symbol",
                    species="mouse",
                    returnall=True
                )

                probe_to_symbol = {}
                for hit in result["out"]:
                    if "symbol" in hit:
                        probe_to_symbol[hit["query"]] = hit["symbol"]

                df_results["Gene"] = df_results["Probe"].map(
                    probe_to_symbol
                ).fillna(df_results["Probe"])
                df_sig["Gene"] = df_sig["Probe"].map(
                    probe_to_symbol
                ).fillna(df_sig["Probe"])

                mapped = df_results["Gene"].ne(df_results["Probe"]).sum()
                st.info(f"{mapped} / {len(df_results)} プローブを遺伝子シンボルに変換しました")

            # 結果表示
            st.markdown("### 解析結果")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("総プローブ数", len(df_results))
            with col2:
                st.metric("有意な変化", len(df_sig))
            with col3:
                st.metric("上昇", len(df_sig[df_sig["Log2FoldChange"] > 0]))

            st.markdown("### 有意に変化したプローブ（上位50件）")
            st.dataframe(
                df_sig[["Gene", "Probe", "Log2FoldChange", "p_value", "p_adj"]].head(50),
                width='stretch'
            )

            st.session_state["df_results"] = df_results
            st.session_state["df_sig"] = df_sig
            st.success("解析完了！Visualization ページに進んでください")

        except Exception as e:
            st.error(f"解析中にエラーが発生しました: {e}")
            st.exception(e)