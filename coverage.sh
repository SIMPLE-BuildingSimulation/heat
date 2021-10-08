rustup component add llvm-tools-preview
cargo install grcov
export RUSTFLAGS="-Zinstrument-coverage"

cargo build

export LLVM_PROFILE_FILE="thermal-%p-%m.profraw"

cargo test

grcov . -s . --binary-path ./target/debug/ -t html --branch --ignore-not-existing -o ./coverage/

rm *.profraw