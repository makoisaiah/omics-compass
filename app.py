# OmicsCompass - メインアプリ
# Streamlit の基本構造について:
# https://docs.streamlit.io/get-started/fundamentals/main-concepts

import streamlit as st

st.set_page_config(
    page_title="OmicsCompass",
    page_icon="🧭",
    layout="wide"
)

st.title("🧭 OmicsCompass")
st.subheader("オミクスデータから経路を探索するツール")

st.markdown("""
### できること（予定）
- 📥 GEO からRNA-seq データを取得
- 🔬 差分発現解析
- 🗺️ パスウェイ解析
- 📊 インタラクティブな可視化
""")

st.info("現在開発中です")
