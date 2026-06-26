# OmicsCompass - 認証モジュール
# Streamlit Secrets を使った認証について:
# https://docs.streamlit.io/develop/concepts/connections/secrets-management
# streamlit-authenticator 0.4.x の使い方:
# https://github.com/mkhorasani/Streamlit-Authenticator

import streamlit as st
import streamlit_authenticator as stauth
import yaml

def load_auth():
    """Streamlit Secrets から認証情報を読み込む"""
    try:
        config = {
            "credentials": {
                "usernames": {
                    k: dict(v) for k, v in st.secrets["credentials"]["usernames"].items()
                }
            },
            "cookie": dict(st.secrets["cookie"])
        }
    except Exception as e:
        st.error(f"認証情報の読み込みに失敗しました: {e}")
        st.info("Streamlit の Secrets に認証情報を設定してください")
        st.stop()

    authenticator = stauth.Authenticate(
        config["credentials"],
        config["cookie"]["name"],
        config["cookie"]["key"],
        config["cookie"]["expiry_days"],
    )
    return authenticator

def login():
    """ログイン画面を表示して認証状態を返す"""
    authenticator = load_auth()

    authenticator.login()

    if st.session_state.get("authentication_status") is False:
        st.error("ユーザー名またはパスワードが違います")
        return False
    elif st.session_state.get("authentication_status") is None:
        st.warning("ユーザー名とパスワードを入力してください")
        return False
    elif st.session_state.get("authentication_status"):
        authenticator.logout(location="sidebar")
        st.sidebar.success(f'ログイン中: {st.session_state["name"]}')
        return True