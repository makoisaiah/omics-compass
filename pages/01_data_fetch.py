# OmicsCompass - GEO データ取得ページ
# NCBI GEO FTP からの補足ファイル取得:
# https://www.ncbi.nlm.nih.gov/geo/info/download.html
# GEOparse の使い方:
# https://geoparse.readthedocs.io/en/latest/

import streamlit as st
import GEOparse
import pandas as pd
import requests
import re
import gzip
import io

st.title("📥 Data Fetch")
st.subheader("GEO からRNA-seq データを取得")

st.markdown("""
GEO（Gene Expression Omnibus）から公開データセットを取得します。
データセット ID（例: GSE21393）を入力してください。
""")

geo_id = st.text_input(
    "GEO データセット ID",
    placeholder="例: GSE190343",
    help="NCBI GEO で検索して見つけた ID を入力してください"
)

if geo_id:
    # GEO ID の形式チェック
    if not re.match(r'^GSE\d+$', geo_id.strip().upper()):
        st.error("GEO ID の形式が正しくありません。GSE から始まる番号を入力してください")
        st.stop()

    geo_id = geo_id.strip().upper()

    # 新しい GEO ID が入力されたらセッションをリセット
    if st.session_state.get("geo_id") != geo_id:
        for key in ["df_expr", "expr_columns", "gse", "df_results", "df_sig", "group_suggestions"]:
            st.session_state.pop(key, None)

    # GEOparse を呼ぶ前に NCBI API で存在確認
    with st.spinner(f"{geo_id} の存在を確認中..."):
        try:
            check = requests.get(
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={"db": "gds", "term": geo_id, "retmode": "json"},
                timeout=10
            )
            result = check.json()
            count = int(result["esearchresult"]["count"])
            if count == 0:
                st.error(f"{geo_id} は GEO に存在しません。ID を確認してください。")
                st.stop()
        except Exception:
            pass  # 確認できなくても続行

    with st.spinner(f"{geo_id} のメタデータを取得中..."):
        try:
            gse = GEOparse.get_GEO(
                geo=geo_id,
                destdir="/tmp",
                silent=True,
                how="brief"
            )
        except Exception as e:
            st.error(f"メタデータの取得に失敗しました: {e}")
            st.stop()

    # 補足ファイルを NCBI FTP から直接探す
    # GEO の FTP 構造: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSExxx nnn/GSExxxx/suppl/
    geo_num = geo_id.replace("GSE", "")
    ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{geo_num[:-3]}nnn/{geo_id}/suppl/"

    with st.spinner("補足ファイルを検索中..."):
        try:
            response = requests.get(ftp_base, timeout=10)
            # HTML からファイル名を抽出
            file_matches = re.findall(
                rf'href="({geo_id}_[^"]+\.(?:txt|csv|tsv)\.gz)"',
                response.text
            )
            supp_files = [ftp_base + f for f in file_matches]
        except Exception:
            supp_files = []

    st.success("メタデータ取得成功！")

    # 基本情報を表示
    st.markdown("### データセット情報")
    col1, col2 = st.columns(2)
    with col1:
        st.metric("サンプル数", len(gse.gsms))
    with col2:
        st.metric("プラットフォーム数", len(gse.gpls))

    st.markdown(f"**タイトル**: {gse.metadata.get('title', ['不明'])[0]}")
    st.markdown(f"**概要**: {gse.metadata.get('summary', ['不明'])[0][:500]}")

    # SuperSeries の場合は SubSeries を案内
    relations = gse.metadata.get("relation", [])
    subseries = [r.split(": ")[-1] for r in relations if "SubSeries of" not in r and "GSE" in r]
    if subseries:
        st.warning(f"このデータセットは SubSeries を含む SuperSeries です。以下の SubSeries を直接指定することをおすすめします: {', '.join(subseries)}")

    # サンプル一覧を表示
    st.markdown("### サンプル一覧")
    sample_info = []
    for gsm_name, gsm in gse.gsms.items():
        sample_info.append({
            "サンプル ID": gsm_name,
            "タイトル": gsm.metadata.get("title", ["不明"])[0],
            "ソース": gsm.metadata.get("source_name_ch1", ["不明"])[0],
        })
    df_samples = pd.DataFrame(sample_info)
    st.dataframe(df_samples, width='stretch')

    # 補足ファイルの確認
    st.markdown("### 補足ファイル")
    
    if supp_files:
        for f in supp_files:
            st.text(f)
        
        # txt.gz または csv.gz ファイルを自動選択
        selected_file = None
        for f in supp_files:
            if any(f.endswith(ext) for ext in [".txt.gz", ".csv.gz", ".tsv.gz"]):
                selected_file = f
                break
        
        if selected_file:
            st.success(f"発現データファイルを検出: {selected_file.split('/')[-1]}")
            
            if st.button("発現データを取得", type="primary"):
                with st.spinner("ダウンロード中..."):
                    try:
                        https_url = selected_file.replace(
                            "ftp://ftp.ncbi.nlm.nih.gov",
                            "https://ftp.ncbi.nlm.nih.gov"
                        )
                        response = requests.get(https_url, timeout=120)
                        
                        with gzip.open(io.BytesIO(response.content), 'rt') as f:
                            df_expr = pd.read_csv(f, sep='\t', index_col=0, low_memory=False)
                        
                        df_expr = df_expr.apply(pd.to_numeric, errors="coerce").dropna()
                        
                        st.success(f"発現データ取得完了: {len(df_expr)} 遺伝子 x {len(df_expr.columns)} サンプル")
                        st.dataframe(df_expr.head(), width='stretch')
                        
                        # セッションに保存
                        st.session_state["gse"] = gse
                        st.session_state["geo_id"] = geo_id
                        st.session_state["df_expr"] = df_expr
                        st.session_state["expr_columns"] = df_expr.columns.tolist()
                        st.success("データを読み込みました。Analysis ページに進んでください。")
                        
                    except Exception as e:
                        st.error(f"発現データの取得に失敗しました: {e}")
        else:
            st.warning("自動検出できる発現データファイルがありません。手動で確認してください。")
    else:
        st.warning("補足ファイルが見つかりませんでした。マイクロアレイデータとして処理します。")
        
        if st.button("マイクロアレイデータを取得", type="primary"):
            with st.spinner("サンプルデータを取得中..."):
                try:
                    # GSM テーブルから発現データを取得
                    expr_data = {}
                    for gsm_name, gsm in gse.gsms.items():
                        if gsm.table is not None and not gsm.table.empty:
                            if "ID_REF" in gsm.table.columns and "VALUE" in gsm.table.columns:
                                table = gsm.table.set_index("ID_REF")["VALUE"]
                                expr_data[gsm_name] = table

                    if expr_data:
                        df_expr = pd.DataFrame(expr_data).dropna()
                        df_expr = df_expr.apply(pd.to_numeric, errors="coerce").dropna()
                        st.success(f"発現データ取得完了: {len(df_expr)} プローブ x {len(df_expr.columns)} サンプル")
                        st.dataframe(df_expr.head(), width='stretch')

                        st.session_state["gse"] = gse
                        st.session_state["geo_id"] = geo_id
                        st.session_state["df_expr"] = df_expr
                        st.success("データを読み込みました。Analysis ページに進んでください。")
                    else:
                        # テーブルが空の場合は gse だけ保存して Analysis に委ねる
                        st.session_state["gse"] = gse
                        st.session_state["geo_id"] = geo_id
                        st.info("サンプルテーブルが空です。Analysis ページで解析を試みてください。")
                except Exception as e:
                    st.error(f"エラー: {e}")