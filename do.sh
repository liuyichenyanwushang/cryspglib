mkdir -p txtfile  # 确保目标目录存在
for file in src/*.rs; do
  cp "$file" "txtfile/$(basename "$file" .rs).txt"
done


for file in src_c/*.c; do
  mv "$file" "src_c/$(basename "$file" .c).txt"
done
