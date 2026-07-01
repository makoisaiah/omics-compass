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

# セッションにデータがあるか確認
if "gse" not in st.session_state:
    st.warning("まず Data Fetch ページでデータを取得してください")
    st.stop()

gse = st.session_state["gse"]
geo_id = st.session_state.get("geo_id", "")

st.markdown(f"**データセット**: {geo_id}")

# サンプル一覧を表示してグループ分けを設定
st.markdown("### グループ分け")
st.markdown("各サンプルをコントロール群・処理群に割り当ててください")

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
    prompt = f"""以下は RNA-seq 実験のサンプル一覧です。
各サンプルがコントロール群か処理群かを判定してください。

サンプル一覧:
{sample_list}

以下の JSON 形式のみで回答してください。他の文章は不要です。
{{
  "サンプルID": "コントロール" または "処理群",
  ...
}}"""

    with st.spinner("LLM が判定中..."):
        try:
            try:
                groq_key = st.secrets["GROQ_API_KEY"]
            except Exception:
                groq_key = None
            result, used_backend = ask_llm(prompt, groq_key)
            st.success(f"{used_backend} で判定しました")

            # JSON をパース
            import re
            json_match = re.search(r'\{.*\}', result, re.DOTALL)
            if json_match:
                suggestions = json.loads(json_match.group())
                st.session_state["group_suggestions"] = suggestions
                # プルダウンの session_state を強制的に更新
                for sample_name, group_value in suggestions.items():
                    st.session_state[f"group_{sample_name}"] = group_value
                st.markdown("**判定結果（確認して必要なら修正してください）**")
                st.json(suggestions)
                st.rerun()
#        else:
#                st.warning("JSON の解析に失敗しました。手動で設定してください。")
        except Exception as e:
            st.error(f"エラー: {e}")

st.divider()

# グループ割り当て UI
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
        # LLM の判定結果があればそれをデフォルトに、なければコントロール
        suggested = suggestions.get(name, "コントロール")
        # セッションにすでに選択済みの値があればそれを優先
        if f"group_{name}" not in st.session_state:
            st.session_state[f"group_{name}"] = suggested
        group = st.selectbox(
            label=name,
            options=["コントロール", "処理群"],
            key=f"group_{name}",
            label_visibility="collapsed"
        )
        group_assignments[name] = group

# 自動判定結果の概要を表示
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
            # まず GSM テーブルから試みる
            expr_data = {}
            for name in sample_names:
                gsm = gse.gsms[name]
                if gsm.table is not None and not gsm.table.empty:
                    table = gsm.table.set_index("ID_REF")["VALUE"]
                    expr_data[name] = table

            if not expr_data:
                # GSM テーブルが空の場合、GSE の補足ファイルから取得を試みる
                # https://geoparse.readthedocs.io/en/latest/
                st.info("個別サンプルのテーブルが空のため、補足ファイルから取得を試みます...")
                
                supp_files = gse.metadata.get("supplementary_file", [])
                
                # SubSeries の補足ファイルも確認
                if not supp_files:
                    import GEOparse
                    for relation in gse.metadata.get("relation", []):
                        if "SubSeries of" not in relation and "SuperSeries of" not in relation:
                            continue
                        # SubSeries の GSE ID を取得
                    # 直接 FTP から取得
                    supp_files = [
                        f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190343/suppl/GSE190343_inhibitor_normalized.txt.gz"
                    ]
                
                ftp_url = None
                for f in supp_files:
                    if f.endswith(".txt.gz") or f.endswith(".csv.gz") or f.endswith(".tsv.gz"):
                        ftp_url = f
                        break
                
                if ftp_url:
                    import io
                    import gzip
                    # FTP URL を HTTPS に変換
                    https_url = ftp_url.replace(
                        "ftp://ftp.ncbi.nlm.nih.gov",
                        "https://ftp.ncbi.nlm.nih.gov"
                    )
                    st.info(f"ダウンロード中: {https_url}")
                    response = requests.get(https_url, timeout=60)
                    
                    with gzip.open(io.BytesIO(response.content), 'rt') as f:
                        df_expr = pd.read_csv(f, sep='\t', index_col=0)
                    
                    # カラム名をサンプル ID に合わせる
                    st.info(f"補足ファイルのカラム: {df_expr.columns.tolist()}")
                    df_expr = df_expr.apply(pd.to_numeric, errors="coerce").dropna()

                    # カラム名（実験名）をサンプル ID に対応付ける
                    # GEO のサンプルタイトルと補足ファイルのカラム名を照合
                    col_to_gsm = {}
                    for gsm_name, gsm in gse.gsms.items():
                        title = gsm.metadata.get("title", [""])[0]
                        for col in df_expr.columns:
                            # タイトルとカラム名の部分一致で対応付け
                            col_clean = col.replace(" ", "_").replace("-", "_").lower()
                            title_clean = title.replace(" ", "_").replace("-", "_").lower()
                            if col_clean in title_clean or title_clean in col_clean:
                                col_to_gsm[col] = gsm_name
                                break

                    st.info(f"カラム対応付け: {col_to_gsm}")

                    if col_to_gsm:
                        df_expr = df_expr.rename(columns=col_to_gsm)
                    else:
                        # 自動対応付けに失敗した場合、順番で対応付け
                        st.warning("自動対応付けに失敗しました。サンプルの順番で対応付けます")
                        col_map = dict(zip(df_expr.columns, sample_names))
                        df_expr = df_expr.rename(columns=col_map)
                        st.info(f"順番による対応付け: {col_map}")
                else:
                    st.error("補足ファイルが見つかりませんでした")
                    st.stop()
            else:
                df_expr = pd.DataFrame(expr_data).dropna()
                df_expr = df_expr.apply(pd.to_numeric, errors="coerce").dropna()

            st.success(f"発現データ取得完了: {len(df_expr)} プローブ")

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

            # 多重検定補正（Benjamini-Hochberg法）
            _, p_adj, _, _ = multipletests(p_values, method="fdr_bh")

            # 結果をデータフレームにまとめる
            df_results = pd.DataFrame({
                "Probe": df_expr.index,
                "Log2FoldChange": log2fc,
                "p_value": p_values,
                "p_adj": p_adj,
            }).sort_values("p_adj")

            # 有意な遺伝子のみ抽出
            df_sig = df_results[
                (df_results["p_adj"] < 0.05) &
                (abs(df_results["Log2FoldChange"]) > 1)
            ].copy()

            # ---- プローブ ID → 遺伝子シンボル変換 ----
            with st.spinner("プローブ ID を遺伝子シンボルに変換中..."):
                mg = mygene.MyGeneInfo()
                probe_ids = df_results["Probe"].tolist()

                result = mg.querymany(
                    probe_ids,
                    scopes="reporter",
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

            # ---- 結果表示 ----
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

            # セッションに結果を保存
            st.session_state["df_results"] = df_results
            st.session_state["df_sig"] = df_sig
            st.success("解析完了！Visualization ページに進んでください")

        except Exception as e:
            st.error(f"解析中にエラーが発生しました: {e}")
            st.exception(e)
