i=12
while [ $i -le 40 ]
do
    [ $i -le 9 ] && ip=0$i || ip=$i
    sed -e "s/XX/"$ip"/" -e "s/YY/"$i"/" cpptraj.in > test.in
    cpptraj -i test.in
    i=`expr $i + 1`
done