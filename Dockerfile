

FROM amazonlinux:2 AS builder

RUN yum -y update && yum clean all && \
    yum -y install \
        wget curl gcc gcc-c++ make openssl openssl-devel perl-core \
        zlib-devel bzip2 bzip2-devel xz xz-devel lzma lzma-devel \
        libcurl libcurl-devel ncurses ncurses-devel \
        which git cmake pkgconfig autoconf automake libtool \
        clang llvm llvm-devel

ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y && \
    rustup default stable

WORKDIR /app

COPY Cargo.toml Cargo.lock ./

COPY src ./src
RUN cargo build --release

FROM amazonlinux:2

WORKDIR /linux_build
COPY --from=builder /app/target/release/isopedia ./isopedia
COPY --from=builder /app/target/release/isopedia-tools ./isopedia-tools

ENTRYPOINT ["./isopedia"]
