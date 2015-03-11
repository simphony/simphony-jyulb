#!/bin/bash
# install JYU-LB
git clone https://github.com/simphony/JYU-LB.git
cd JYU-LB
make
ln -s /usr/local/bin/jyu_lb_isothermal3D.exe bin/jyu_lb_isothermal3D.exe
