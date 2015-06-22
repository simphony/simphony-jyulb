#!/bin/bash
# install JYU-LB
prefix=/usr/local/bin/
args=$(getopt -l "prefix:" -o "s:h" -- "$@")

eval set -- "$args"

while [ $# -ge 1 ]; do
        case "$1" in
                --)
                    # No more options left.
                    shift
                    break
                   ;;
                -p|--prefix)
                        prefix="$2"
                        shift
                        ;;
                -h)
                        echo "Usage: ./install_external.sh [--prefix=somewhere]"
                        exit 0
                        ;;
        esac

        shift
done
echo "Installing to $prefix"

# Make a temporary directory for checkout
temp_checkout_dir=`mktemp -d`

git clone https://github.com/simphony/JYU-LB.git $temp_checkout_dir
cd $temp_checkout_dir
make -j 2

echo "Moving the executable to $prefix/jyu_lb_isothermal.exe"
# Move the executable to the desired path
mv  $(pwd)/bin/jyu_lb_isothermal.exe "$prefix/jyu_lb_isothermal.exe"
cd ..

echo "Cleaning up.."
# Delete the temprary directory
rm -rf $temp_checkout_dir
echo "Bye"
