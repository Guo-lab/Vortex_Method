make clean

make vortex

cd vtk

expect <<EOF
spawn ../vortexMethod/vortex2D.exe

expect "input log_2(number of grid points)"
send "6\r"

expect "input test = 1,2, other"
send "3\r"

expect "input particle refinement factor"
send "8\r"

expect "enter stopping time"
send "10\r"

expect eof
EOF
