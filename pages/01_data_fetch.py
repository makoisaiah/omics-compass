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
st.subheader("GEO からRNA-seq・マイクロアレイデータを取得")

geo_id = st.text_input(
    "GEO データセット ID",
    placeholder="例: GSE21393",
    help="NCBI GEO で検索して見つけた ID を入力してください"
)

if geo_id:
    if not re.match(r'^GSE\d+$', geo_id.strip().upper()):
        st.error("GSE から始まる番号を入力してください（例: GSE21393）")
        st.stop()

    geo_id = geo_id.strip().upper()

    # 新しい GEO ID が入力されたらセッションをリセット
    if st.session_state.get("geo_id") != geo_id:
        for key in ["df_expr", "gse", "df_results", "df_sig", "group_suggestions"]:
            st.session_state.pop(key, None)

    # NCBI API で存在確認
    with st.spinner(f"{geo_id} を確認中..."):
        try:
            check = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={"db": "gds", "term": geo_id, "retmode": "json"},
                timeout=10
            )
            count = int(check.json()["esearchresult"]["count"])
            if count == 0:
                st.error(f"{geo_id} は GEO に存在しません")
                st.stop()
        except Exception:
            pass

    # メタデータ取得
    with st.spinner(f"{geo_id} のメタデータを取得中..."):
        try:
            gse = GEOparse.get_GEO(
                geo=geo_id,
                destdir="/tmp",
                silent=True,
                include_data=False
            )
        except Exception as e:
            st.error(f"メタデータの取得に失敗しました: {e}")
            st.stop()

    # 基本情報を表示
    st.success("メタデータ取得成功！")
    st.markdown(f"**タイトル**: {gse.metadata.get('title', ['不明'])[0]}")
    st.markdown(f"**概要**: {gse.metadata.get('summary', ['不明'])[0][:300]}")

    col1, col2 = st.columns(2)
    with col1:
        st.metric("サンプル数", len(gse.gsms))
    with col2:
        st.metric("プラットフォーム数", len(gse.gpls))

    # SuperSeries チェック
    relations = gse.metadata.get("relation", [])
    subseries = [
        r.split("SuperSeries of: ")[-1].strip()
        for r in relations
        if "SuperSeries of:" in r
    ]
    if subseries:
        st.warning(f"このデータセットは SuperSeries です。以下の SubSeries を直接指定することをおすすめします。")
        for s in subseries:
            st.code(s)
        st.stop()

    # サンプル一覧
    st.markdown("### サンプル一覧")
    sample_info = []
    for gsm_name, gsm in gse.gsms.items():
        sample_info.append({
            "サンプル ID": gsm_name,
            "タイトル": gsm.metadata.get("title", ["不明"])[0],
            "ソース": gsm.metadata.get("source_name_ch1", ["不明"])[0],
        })
    st.dataframe(pd.DataFrame(sample_info), width='stretch')

    # データ種別を判定
    series_type = gse.metadata.get("type", [""])[0].lower()
    is_rnaseq = "high throughput sequencing" in series_type

    # 補足ファイルを FTP から探す
    geo_num = geo_id.replace("GSE", "")
    ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{geo_num[:-3]}nnn/{geo_id}/suppl/"
    supp_files = []
    try:
        resp = requests.get(ftp_base, timeout=10)
        matches = re.findall(
            rf'href="({geo_id}_[^"]+\.(?:txt|csv|tsv)\.gz)"',
            resp.text
        )
        supp_files = [ftp_base + f for f in matches]
    except Exception:
        pass

    # NCBI 自動生成カウント行列を探す（FTP に補足ファイルがない場合の代替）
    # https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html
    if not supp_files:
        ncbi_genomes = [
            "GRCh38.p13_NCBI",
            "GRCh38.p14_NCBI", 
            "GRCm38.p6_NCBI",
            "GRCm39_NCBI",
        ]
        for genome in ncbi_genomes:
            ncbi_url = (
                f"https://www.ncbi.nlm.nih.gov/geo/download/"
                f"?type=rnaseq_counts&acc={geo_id}&format=file"
                f"&file={geo_id}_raw_counts_{genome}.tsv.gz"
            )
            try:
                check = requests.head(ncbi_url, timeout=5)
                if check.status_code == 200:
                    supp_files = [ncbi_url]
                    st.info(f"NCBI 自動生成カウント行列を検出: {genome}")
                    break
            except Exception:
                pass

    st.markdown("### データ取得")

    if supp_files:
        # RNA-seq: 補足ファイルから取得
        st.info(f"補足ファイルを検出しました（RNA-seq）")
        for f in supp_files:
            st.text(f.split("/")[-1])

        selected_file = supp_files[0]

        if st.button("発現データを取得", type="primary"):
            with st.spinner("ダウンロード中..."):
                try:
                    response = requests.get(selected_file, timeout=120)
                    with gzip.open(io.BytesIO(response.content), 'rt') as f:
                        df_expr = pd.read_csv(f, sep='\t', index_col=0, low_memory=False, comment='#')
                    with gzip.open(io.BytesIO(response.content), 'rt') as f:
                        df_expr = pd.read_csv(f, sep='\t', index_col=0, low_memory=False)
                    df_expr = df_expr.apply(pd.to_numeric, errors="coerce").dropna()

                    st.success(f"取得完了: {len(df_expr)} 遺伝子 x {len(df_expr.columns)} サンプル")
                    st.dataframe(df_expr.head(), width='stretch')

                    st.session_state["gse"] = gse
                    st.session_state["geo_id"] = geo_id
                    st.session_state["df_expr"] = df_expr
                    st.success("Analysis ページに進んでください。")
                except Exception as e:
                    st.error(f"取得に失敗しました: {e}")

    else:
        # マイクロアレイ: GSM テーブルから取得
        st.info("補足ファイルなし。マイクロアレイデータとして取得します。")

        if st.button("マイクロアレイデータを取得", type="primary"):
            with st.spinner("サンプルデータを取得中（時間がかかる場合があります）..."):
                try:
                    # include_data=False で取得したので改めてフルで取得
                    gse_full = GEOparse.get_GEO(
                        geo=geo_id,
                        destdir="/tmp",
                        silent=True
                    )
                    expr_data = {}
                    for gsm_name, gsm in gse_full.gsms.items():
                        if gsm.table is not None and not gsm.table.empty:
                            if "ID_REF" in gsm.table.columns and "VALUE" in gsm.table.columns:
                                table = gsm.table.set_index("ID_REF")["VALUE"]
                                expr_data[gsm_name] = table

                    if expr_data:
                        df_expr = pd.DataFrame(expr_data).dropna()
                        df_expr = df_expr.apply(pd.to_numeric, errors="coerce").dropna()

                        st.success(f"取得完了: {len(df_expr)} プローブ x {len(df_expr.columns)} サンプル")
                        st.dataframe(df_expr.head(), width='stretch')

                        st.session_state["gse"] = gse_full
                        st.session_state["geo_id"] = geo_id
                        st.session_state["df_expr"] = df_expr
                        st.success("Analysis ページに進んでください。")
                    else:
                        st.error("発現データが見つかりませんでした")
                except Exception as e:
                    st.error(f"取得に失敗しました: {e}")