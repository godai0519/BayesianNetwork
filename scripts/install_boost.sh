#!/bin/bash

readonly tmp_dir=/tmp
readonly work_dir=${tmp_dir}/boost_${BOOST_VERSION//[.]/_}

function download_source() {
    local local_path=${tmp_dir}/boost.tar.gz
    local url_prefix=http://sourceforge.net/projects/boost/files/boost
    local url=$url_prefix/${BOOST_VERSION}/boost_${BOOST_VERSION//[.]/_}.tar.gz/download
    wget -q -O $local_path $url || return 1
    tar xzf $local_path -C ${tmp_dir} || return 1
    return 0
}

function install_boost {
    (
        cd ${work_dir} && \
        ./bootstrap.sh --prefix=/usr --libdir=/usr/lib && \
        ./b2  --with-test --stagedir=stage --prefix=/usr --libdir=/usr/lib -j2
    )
    return $?
}

download_source || exit 1
install_boost || exit 1
exit 0
