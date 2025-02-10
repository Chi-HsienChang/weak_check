# weak_check
weak_check

chmod +x run_no_weak.sh
./run_no_weak.sh


g++ -o a.out a.cpp -fopenmp
./a.out 4 emu 125 > a.txt
./b.out > b.txt

check storage:
du -sh .


./a_S_L.out 4 emu 25 > ell_4_S_L.txt

python3 read.py ell_4_S_L.txt

sed -n '1,100p' ell_4_S_L.txt 