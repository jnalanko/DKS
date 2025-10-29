set -xue

cargo build --release --target x86_64-unknown-linux-musl

VERSION=$(cargo metadata --no-deps --format-version=1 | jq -r ".packages[0].version")
TAG="v$VERSION"

RELEASE_DIR=dks-${TAG}-linux-x86_64
mkdir -p $RELEASE_DIR
cp README.md LICENSE target/x86_64-unknown-linux-musl/release/dks $RELEASE_DIR
tar -cvzf $RELEASE_DIR.tar.gz $RELEASE_DIR


echo git tag -a "$TAG" -m "DKS $TAG"
echo git push origin "$TAG"

