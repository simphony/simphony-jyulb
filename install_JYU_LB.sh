#!/bin/bash
# install JYU-LB modeling engine
git clone https://github.com/simphony/JYU-LB.git
cd JYU-LB
make
cd ..
mv JYU-LB/bin/jyu_lb_isothermal3D.exe .
