for file in *TCcounts; do awk '{sum+=$2}END{print sum}' $file > $file.total; done
