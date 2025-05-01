make clean

make vortex

cd vtk

expect <<EOF
spawn ../vortexMethod/vortex2D.exe

expect "input log_2(number of grid points)"
send "7\r"

expect "input test = 1,2, other"
send "5\r"

expect "input particle refinement factor"
send "2\r"

expect "enter stopping time"
send "12.5\r"

expect eof
EOF
