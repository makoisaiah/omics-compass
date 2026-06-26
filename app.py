# OmicsCompass - メインアプリ
# Streamlit のマルチページアプリについて:
# https://docs.streamlit.io/develop/concepts/multipage-apps/overview

import streamlit as st

st.set_page_config(
    page_title="OmicsCompass",
    page_icon="🧭",
    layout="wide"
)

st.title("🧭 OmicsCompass")
st.subheader("オミクスデータから経路を探索するツール")

st.markdown("""
### ナビゲーション
左のサイドバーからページを選んでください。

| ページ | 内容 |
|--------|------|
| 📥 Data Fetch | GEO からRNA-seq データを取得 |
| 🔬 Analysis | 差分発現解析 |
| 📊 Visualization | パスウェイ解析と可視化 |
""")

st.info("現在開発中です")