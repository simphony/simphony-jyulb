#!/bin/bash
# install JYU-LB
git clone https://github.com/simphony/JYU-LB.git
cd JYU-LB
make
ln -s  $(pwd)/bin/jyu_lb_isothermal3D.exe /usr/local/bin/jyu_lb_isothermal3D.exe
