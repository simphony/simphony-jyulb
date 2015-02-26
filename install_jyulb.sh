#!/bin/bash
# install JYU-LB
git clone https://github.com/simphony/JYU-LB.git
cd JYU-LB
make
cd ..
mv JYU-LB/bin/jyu_lb_isothermal3D.exe .
