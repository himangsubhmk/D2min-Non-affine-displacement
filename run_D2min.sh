for i in {1..99}
do
    j=$((i-1))

    ./a.out 10000 "$j" "$i" 
done
