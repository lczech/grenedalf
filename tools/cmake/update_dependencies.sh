#!/bin/bash

# grenedalf - Genome Analyses of Differential Allele Frequencies
# Copyright (C) 2020-2022 Lucas Czech
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Lucas Czech <lczech@carnegiescience.edu>
# Department of Plant Biology, Carnegie Institution For Science
# 260 Panama Street, Stanford, CA 94305, USA

####################################################################################################
#    Update Git Submodule Commit Hashes
####################################################################################################

# Here, we update the commit hashes of the download scripts for our dependencies.
# The script uses git to get the hases that the submodules currently point to,
# and replaces them in the cmake download scripts.

# Change to main project dir. This ensures that the script can be called from any directory.
cd "$(dirname "$0")/../.."

function get_commit_hash() {
    local libname=${1}

    # Disect the git submodule output to get the commit that is currently checked out.
    # This output is of the form ".<hash> libs/genesis (v0.19.0-14-g5b77dba)",
    # where the first dot is a space, plus sign or some other git status symbol.
    # This also sets the "return" to the outside scope... bash...
    commit_hash=`git submodule status | grep ${libname} | sed "s/.\([0-9a-f]*\) .*/\1/g" | tr -d "\n"`

    # Use two methods for checking.
    # This fails if the submodule is checked out, but not commited.
    # In that case, we only print a warning.
    local branch=`git rev-parse --abbrev-ref HEAD`
    local hash_a=`git ls-files -s libs/${libname} | cut -d" " -f 2`
    local hash_b=`git ls-tree ${branch} libs/${libname} | awk -F "[ \t]" '{print $3}'`

    if [[ ${hash_a} != ${hash_b} ]]; then
        echo -e "\e[31mProblem with commit hash for ${libname}: ${hash_a} != ${hash_b}\e[0m"
    fi
    if [[ ${hash_a} != ${commit_hash} ]]; then
        echo -e "\e[31mNot yet commited submodule ${libname}: ${commit_hash} > ${hash_a}\e[0m"
    fi
}

function update_commit_hash() {
    local libname=${1}

    get_commit_hash ${libname}
    echo "${libname} @ ${commit_hash}"

    LINE="SET( ${libname}_COMMIT_HASH \"${commit_hash}\" ) #${libname}_COMMIT_HASH#"
    sed -i "s/.*#${libname}_COMMIT_HASH#/${LINE}/g" CMakeLists.txt
}

####################################################################################################
#    Submodule Dependencies
####################################################################################################

update_commit_hash "CLI11"
update_commit_hash "genesis"
