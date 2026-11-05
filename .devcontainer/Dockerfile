# vormorph の開発環境（GitHub Codespaces / VS Code Dev Containers）
#
# - Rust：ベースイメージの rustup を使う。ツールチェーンは rust-toolchain.toml に従う
# - Python：uv で tools/ の環境を作る（Python 本体も uv が用意する）
# - gmsh（PyPI）：OpenGL と X11 の共有ライブラリを必要とする
FROM mcr.microsoft.com/devcontainers/rust:1-trixie

# gmsh の libgmsh.so が必要とする共有ライブラリ。
# CI の ubuntu-latest には最初から入っているものも、コンテナでは明示的に入れる。
# tmux は長い計算をターミナルから切り離して続けるために使う。
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y --no-install-recommends \
        libglu1-mesa \
        libxft2 \
        libxcursor1 \
        libxinerama1 \
        libgomp1 \
        tmux \
    && rm -rf /var/lib/apt/lists/*

# uv は公式のイメージから実行ファイルだけを取り出す（版を固定する）
COPY --from=ghcr.io/astral-sh/uv:0.13.7 /uv /uvx /usr/local/bin/
