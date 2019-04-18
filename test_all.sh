#!/usr/bin/env sh
echo_green(){
    echo -e "\e[1;42m$1\e[0m"
}

set -e




cd Task4
make test

 echo_green 'All tests done!'