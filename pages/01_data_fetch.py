# OmicsCompass - GEO データ取得ページ
# GEOparse の使い方:
# https://geoparse.readthedocs.io/en/latest/
# NCBI GEO データベース:
# https://www.ncbi.nlm.nih.gov/geo/

import streamlit as st
import GEOparse
import pandas as pd

st.title("📥 Data Fetch")
st.subheader("GEO からRNA-seq データを取得")

st.markdown("""
GEO（Gene Expression Omnibus）から公開データセットを取得します。
データセット ID（例: GSE123456）を入力してください。
""")

# データセット ID の入力
geo_id = st.text_input(
    "GEO データセット ID",
    placeholder="例: GSE190346",
    help="NCBI GEO で検索して見つけた ID を入力してください"
)

if geo_id:
    with st.spinner(f"{geo_id} を取得中..."):
        try:
            # GEO からデータを取得
            # ダウンロードしたファイルは /tmp に一時保存（ローカルに残さない）
            gse = GEOparse.get_GEO(
                geo=geo_id,
                destdir="/tmp",
                silent=True
            )

            st.success(f"取得成功！")

            # データセットの基本情報を表示
            st.markdown("### データセット情報")
            col1, col2 = st.columns(2)
            with col1:
                st.metric("サンプル数", len(gse.gsms))
            with col2:
                st.metric("プラットフォーム数", len(gse.gpls))

            st.markdown(f"**タイトル**: {gse.metadata.get('title', ['不明'])[0]}")
            st.markdown(f"**概要**: {gse.metadata.get('summary', ['不明'])[0][:500]}...")

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

            # セッションにデータを保存（次のページで使えるように）
            st.session_state["gse"] = gse
            st.session_state["geo_id"] = geo_id
            st.success("データを読み込みました。Analysis ページに進んでください。")

        except Exception as e:
            st.error(f"取得に失敗しました: {e}")
            st.info("GEO ID が正しいか確認してください。例: GSE190346")