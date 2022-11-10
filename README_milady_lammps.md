
```
pluma ~/README_milady_lammps.md
```

```
sudo bash
yum install virtualenv
```

### git push lammps stable modifications for milady

```
pwd
    /home/catA/wambeke/MLD_LAMMPS/ml_cv1/lammps
git remote -v
    origin	/home/catA/wambeke/MLD_LAMMPS/lammps_stable (fetch)
    origin	/home/catA/wambeke/MLD_LAMMPS/lammps_stable (push)
cd /home/catA/wambeke/MLD_LAMMPS
git clone git@github.com:wambeke/lammps.git lammps_stable_cv
cd lammps_stable_cv
git checkout stable
git checkout -b stable_milady_cv
meld ../lammps_stable_cv ../ml_cv1/lammps



```

### compilation doc

```
https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork

Configuring a remote for a fork
You must configure a remote that points to the upstream repository in Git 
to sync changes you make in a fork with the original repository. 
This also allows you to sync changes made in the original repository with the fork.

$ git remote -v
> origin  https://github.com/YOUR_USERNAME/YOUR_FORK.git (fetch)
> origin  https://github.com/YOUR_USERNAME/YOUR_FORK.git (push)

Specify a new remote upstream repository that will be synced with the fork.

$ git remote add upstream https://github.com/ORIGINAL_OWNER/ORIGINAL_REPOSITORY.git

Verify the new upstream repository you've specified for your fork.

$ git remote -v
> origin    https://github.com/YOUR_USERNAME/YOUR_FORK.git (fetch)
> origin    https://github.com/YOUR_USERNAME/YOUR_FORK.git (push)
> upstream  https://github.com/ORIGINAL_OWNER/ORIGINAL_REPOSITORY.git (fetch)
> upstream  https://github.com/ORIGINAL_OWNER/ORIGINAL_REPOSITORY.git (push)
```

```
# https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork

Syncing a fork branch from the command line
Before you can sync your fork with an upstream repository, 
you must configure a remote that points to the upstream repository in Git.

    Open Terminal.
    Change the current working directory to your local project.
    Fetch the branches and their respective commits from the upstream repository. 
    Commits to BRANCHNAME will be stored in the local branch upstream/BRANCHNAME.

    $ git fetch upstream
    > remote: Counting objects: 75, done.
    > remote: Compressing objects: 100% (53/53), done.
    > remote: Total 62 (delta 27), reused 44 (delta 9)
    > Unpacking objects: 100% (62/62), done.
    > From https://github.com/ORIGINAL_OWNER/ORIGINAL_REPOSITORY
    >  * [new branch]      main     -> upstream/main

    Check out your fork's local default branch - in this case, we use main.

    $ git checkout main
    > Switched to branch 'main'

    Merge the changes from the upstream default branch 
    - in this case, upstream/main - into your local default branch. 
    This brings your fork's default branch into sync with the upstream repository,
    without losing your local changes.

    $ git merge upstream/main
    > Updating a422352..5fdff0f
    > Fast-forward
    >  README                    |    9 -------
    >  README.md                 |    7 ++++++
    >  2 files changed, 7 insertions(+), 9 deletions(-)
    >  delete mode 100644 README
    >  create mode 100644 README.md

    If your local branch didn't have any unique commits, 
    Git will perform a fast-forward. 
    For more information, see Basic Branching and Merging in the Git documentation.

    $ git merge upstream/main
    > Updating 34e91da..16c56ad
    > Fast-forward
    >  README.md                 |    5 +++--
    >  1 file changed, 3 insertions(+), 2 deletions(-)

    If your local branch had unique commits, you may need to resolve conflicts. 
    For more information, see "Addressing merge conflicts."
```


```
export MLD_LAMMPS_ROODIR=~/MLD_LAMMPS
mkdir ${MLD_LAMMPS_ROODIR}
cd ${MLD_LAMMPS_ROODIR}

# get original release lammps
export MLD_LAMMPS_ORI_SRCDIR=${MLD_LAMMPS_ROODIR}/lammps
git clone --branch release https://github.com/lammps/lammps.git  # as ${MLD_LAMMPS_ORI_SRCDIR}

export MLD_LAMMPS_SRCDIR=${MLD_LAMMPS_ROODIR}/milady_lammps
git clone --branch 2022 git@github.com:mcmarinica/milady_lammps.git

cd milady_lammps/
git branch -a

git remote add upstream git@github.com:lammps/lammps.git
git remote set-url --push upstream DISABLE
git remote -v
git fetch upstream

# https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork
# needs a git checkout and a git merge
# git checkout main
# git merge upstream/main

git st


cd ${MLD_LAMMPS_ROODIR}
git clone git@github.com:mcmarinica/milady_lammps.git milady_lammps_ori

meld milady_lammps* &

cd ${MLD_LAMMPS_ROODIR}
cd milady_lammps
find . -name "*pdf" -type f
cd doc
make html
firefox html/index.html &
make pdf
evince Manual.pdf &
```

### compilation make

OK

```
export MLD_LAMMPS_SRCDIR=${MLD_LAMMPS_ROODIR}/milady_lammps


cd ${MLD_LAMMPS_SRCDIR}/src  # change to main LAMMPS source folder

make -j 10 serial  # build a serial LAMMPS executable cd using GNU g++

module avail
  mp-x86_64  mpi/mpich-x86_64  mpi/openmpi-x86_64  

module load mpi/openmpi-x86_64
make -j 10 mpi     # build a parallel LAMMPS executable with MPI

make               # see a variety of make options
```


### compilation cmake

OK if cmake/etc present as merge should do, or you do it a mano.

```
cd ${MLD_LAMMPS_SRCDIR}  # change to main LAMMPS source folder
rm -rf build
mkdir build
cd build
cmake ../cmake
make -C /home/catA/wambeke/MLD_LAMMPS/milady_lammps/src purge
cd ..
rm -rf build
mkdir build
cd build
cmake ../cmake
cmake --build .
make install
```

Problem if not 

```

cd ${MLD_LAMMPS_SRCDIR}  # change to main LAMMPS source folder

cmake --version
rm -rf build
mkdir build
cd build
cmake ../cmake # -> works only if copy ${MLD_LAMMPS_ROODIR}/lammps/cmake/etc to ${MLD_LAMMPS_SRCDIR}/cmake
cmake --build .   


cmake ../cmake # -> do not works
      cmake ../cmake
      -- The CXX compiler identification is GNU 11.3.1
...
      -- Generating package headers...
      -- Generating lmpinstalledpkgs.h...
      CMake Error: File /home/catA/wambeke/MLD_LAMMPS/milady_lammps/cmake/etc/profile.d/lammps.sh.in does not exist.
      CMake Error at CMakeLists.txt:734 (configure_file):
        configure_file Problem configuring file


      CMake Error: File /home/catA/wambeke/MLD_LAMMPS/milady_lammps/cmake/etc/profile.d/lammps.csh.in does not exist.
      CMake Error at CMakeLists.txt:735 (configure_file):
        configure_file Problem configuring file


...
    See also "/home/catA/wambeke/MLD_LAMMPS/milady_lammps/build/CMakeFiles/CMakeError.log".
```

### set lammps in ml cv

Create https://github.com/wambeke/ml as fork of https://github.com/mcmarinica/ml
in github ihm


```
cd ${MLD_LAMMPS_ROODIR}
export MLD_ML_CV_SRCDIR=${MLD_LAMMPS_ROODIR}/ml_cv
git clone --branch master_cv git@github.com:wambeke/ml.git ml_cv  # wambeke is fork of mcmarinica/ml
cd ${MLD_ML_CV_SRCDIR}

# add lammps github
git submodule add git@github.com:lammps/lammps.git lammps
    Clonage dans '/export/home/catA/wambeke/MLD_LAMMPS/ml_cv/lammps'...
    remote: Enumerating objects: 338148, done.
    remote: Counting objects: 100% (78/78), done.
    remote: Compressing objects: 100% (53/53), done.
    remote: Total 338148 (delta 37), reused 38 (delta 25), pack-reused 338070
    Réception d'objets: 100% (338148/338148), 616.10 Mio | 14.64 Mio/s, fait.
    Résolution des deltas: 100% (285099/285099), fait.

git st
    Sur la branche master_cv
    Votre branche est à jour avec 'origin/master_cv'.
    Modifications qui seront validées :
      (utilisez "git restore --staged <fichier>..." pour désindexer)
        nouveau fichier : .gitmodules
        nouveau fichier : lammps
        
cd lammps
@is246206.../lammps(develop)>git branch -a
  * develop
    remotes/origin/HEAD -> origin/develop
    remotes/origin/amoeba-gpu
    remotes/origin/coulomb-refactoring
    remotes/origin/develop
    remotes/origin/distributed-grids
    remotes/origin/nwchem
    remotes/origin/release
    remotes/origin/stable
    remotes/origin/subversion

git checkout release

git add -A
@is246206.../ml_cv(master_cv)>git st
    Sur la branche master_cv
    Votre branche est à jour avec 'origin/master_cv'.

    Modifications qui seront validées :
      (utilisez "git restore --staged <fichier>..." pour désindexer)
        nouveau fichier : .gitmodules
        nouveau fichier : lammps

@is246206.../ml_cv(master_cv)>git commit -m "submodule add git@github.com:lammps/lammps.git"
    [master_cv 01ce637] submodule add git@github.com:lammps/lammps.git
     2 files changed, 4 insertions(+)
     create mode 100644 .gitmodules
     create mode 160000 lammps

git push
    Énumération des objets: 4, fait.
    Décompte des objets: 100% (4/4), fait.
    Compression par delta en utilisant jusqu'à 20 fils d'exécution
    Compression des objets: 100% (3/3), fait.
    Écriture des objets: 100% (3/3), 409 octets | 409.00 Kio/s, fait.
    Total 3 (delta 1), réutilisés 0 (delta 0), réutilisés du pack 0
    remote: Resolving deltas: 100% (1/1), completed with 1 local object.
    To github.com:wambeke/ml.git
       26760f2..01ce637  master_cv -> master_cv

```

Verif 

```
cd ${MLD_LAMMPS_ROODIR}
git clone --branch master_cv git@github.com:wambeke/ml.git ml_cv_TMP
  -> ml_cv_TMP/lammps VIDE
git submodule update --init --recursive
  -> ml_cv_TMP/lammps non vide

git clone --branch master_cv --recurse-submodules git@github.com:wambeke/ml.git ml_cv_TMP2
  -> ml_cv_TMP/lammps non vide 
    * e745c8aac4 - (HEAD, tag: patch_15Sep2022, origin/release) Merge pull request #3446 from akohlmey/next-patch-release (il y a 6 semaines) <Axel Kohlmeyer>
    (HEAD détachée sur e745c8aac4)
  
git clone --branch master_cv --recurse-submodules git@github.com:wambeke/ml.git ml_cv_TMP3
```

### DIFFERENCES

```
cd /home/catA/wambeke/MLD_LAMMPS
ls
  lammps_develop
  lammps_release
  lammps_stable
  milady_lammps

export ii="lammps_develop lammps_release lammps_stable milady_lammps_ori"
for i in $ii; do
  echo $i
  tree -d $i > tree_$i.tmp
done

meld tree_lammps_*.tmp &
meld tree_lammps_stable.tmp tree_milady_lammps_ori.tmp &
meld lammps_stable milady_lammps_ori &
meld 


# examples/PACKAGES/milady src/ML-MILADY src/USER-MILADY lib/milady
export SRC=milady_lammps_ori
export TRG=lammps_stable
cp -r $SRC/examples/PACKAGES/milady $TRG/examples/PACKAGES/.
cp -r $SRC/src/ML-MILADY $TRG/src/.
cp -r $SRC/src/USER-MILADY $TRG/src/.
cp -r $SRC/lib/milady $TRG/lib/.
cp $SRC/README.md $TRG/.
cp $SRC/lammps_options.cmake $TRG/.

meld $SRC/src/version.h $TRG/src/version.h &

for i in ML-MILADY USER-MILADY; do
  find $SRC -type f -name "*" -exec fgrep -Hni $i {} \;
  find $TRG -type f -name "*" -exec fgrep -Hni $i {} \;
done

fgrep ML-MILADY $SRC/src/*
```

# AS milady_lammps_ori/README.md

## MILADY PRIVATE FORK

Directories differing from main LAMMPS repo:
```bash
src/ML-MILADY/
lib/milady/
examples/PACKAGES/milady
```

## Keeping up-to-date with LAMMPS: private fork hack

You should run the following git commands to fetch from the public LAMMPS branch:
```bash
git remote add upstream git@github.com:lammps/lammps.git
git remote set-url --push upstream DISABLE
```

Then run
```bash
git remote -v
```
You should see
```bash
origin	git@github.com:tomswinburne/milady_lammps.git (fetch)
origin	git@github.com:tomswinburne/milady_lammps.git (push)
upstream	git@github.com:lammps/lammps.git (fetch)
upstream	DISABLE (push)
```

To update (WARNING- this may require a large download!!)

```bash
git fetch upstream
```


## Prepare ml_cv

See https://github.com/mcmarinica/milady_lammps/README.md

```
export MLD_LAMMPS_ROODIR=~/MLD_LAMMPS
export SRC0=ml_cv
export SRC1=milady_lammps_ori
export SRC2=lammps_stable
export TRG1=ml_cv1

cd $MLD_LAMMPS_ROODIR
rm -rf $TRG1

# first option recurse-submodules
# git clone --branch master_cv --recurse-submodules ml_cv $TRG1

# second option git clone stable
git clone --branch master_cv ml_cv $TRG1
cd $TRG1
rm -rf lammps
# git clone --branch stable git@github.com:lammps/lammps.git
git clone --branch stable $MLD_LAMMPS_ROODIR/lammps_stable lammps

cd lammps
git his1

cd $MLD_LAMMPS_ROODIR

cp -r $SRC1/examples/PACKAGES/milady $TRG1/lammps/examples/PACKAGES/.
cp -r $SRC1/src/ML-MILADY $TRG1/lammps/src/.
cp -r $SRC1/src/USER-MILADY $TRG1/lammps/src/.
cp -r $SRC1/lib/milady $TRG1/lammps/lib/.
cp $SRC1/README.md $TRG1/lammps/.
cp $SRC1/lammps_options.cmake $TRG1/lammps/.

cd $MLD_LAMMPS_ROODIR/$TRG1
git st

cd $MLD_LAMMPS_ROODIR/$TRG1/lammps
git st

```


### Compilation milady_lammps_cmake


```
cd /export/home/catA/wambeke/MLD_LAMMPS
cp -rf milady_lammps milady_lammps_cmake
# cd milady_lammps_cmake
# git his1

export TRG7=milady_lammps_cmake
cd $MLD_LAMMPS_ROODIR/$TRG7


rm -rf build
mkdir build
cd build

bash

# on is246206 compile milady and run ctest and run python tests
# rename and modification and use at your own risk

# bash
# export PATH=/export/home/catA/wambeke/MLD/MILADY/scripts:${PATH}
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_intel_is246206.bash
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_intel_is246206.bash | grep -e '-- ..' | grep MLD_MKL_LIB
# f_stuff_intel_is246206


cmake --version
    cmake version 3.20.5

export LOC_ROODIR="/export/home/catA/wambeke"
export MLD_OAP_ROODIR="/export/home/catA/intel"
export MLD_OAP_INSDIR="${MLD_OAP_ROODIR}/oneapi"
  
export CONFIG_SETVARS=../config_setvars_oneapi.tmp
cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
mpi=latest
itac=latest
EOT

cat ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
mpi=latest
itac=latest

source ${MLD_OAP_INSDIR}/setvars.sh --config=${CONFIG_SETVARS}
 
    :: initializing oneAPI environment ...
       bash: BASH_VERSION = 5.1.0(1)-release
       args: Using "$@" for setvars.sh arguments: --config=../config_setvars_oneapi.tmp
    :: compiler -- latest
    :: mkl -- latest
    :: mpi -- latest
    :: itac -- latest
    :: oneAPI environment initialized ::

CMAKE_C_COMPILER C compiler to be used for compilation (default: system specific, gcc on Linux)
CMAKE_CXX_COMPILER C++ compiler to be used for compilation (default: system specific, g++ on Linux)
CMAKE_Fortran_COMPILER Fortran compiler to be used for compilation (default: system specific, gfortran on Linux)

export CCOMP=$(which icc)
export CXXCOMP=$(which icpc)
export FORTRANCOMP=$(which ifort)
  #export FC=$(which ifort)
  #export CC=$(which icc)
  
export MKL_ROOT=${MKLROOT}

icc -help | head
icpc -help | head
ifort -help | head

export TRG6=ml_cv1/lammps

cd $MLD_LAMMPS_ROODIR/$TRG7
rm -rf build ; mkdir build ; cd build

cat ../cmake/presets/basic.cmake # cvw version stable
    # preset that turns on just a few, frequently used packages
    # this will be compiled quickly and handle a lot of common inputs.

    set(ALL_PACKAGES KSPACE MANYBODY MOLECULE RIGID)

    foreach(PKG ${ALL_PACKAGES})
      set(PKG_${PKG} ON CACHE BOOL "" FORCE)
    endforeach()
    
cat ../lammps_options.cmake  # marinica version
    # FOR LAMMPS CMAKE BUILD

    # set installation location
    set(CMAKE_INSTALL_PREFIX "$ENV{PREFIX}")

    # enforce c++11 standards
    set(CCFLAGS -g -O3 -std=c++11)

    # compile a binary and shared library
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

    # allow error messages (very useful)
    set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)

    # minimal packages to run example (MANYBODY, EXTRA-FIX) and
    # generate new pathway (REPLICA for "fix neb")
    # Also include some "ML" potentials (PACE and SNAP)
    set(ALL_PACKAGES MANYBODY EXTRA-FIX REPLICA ML-MILADY ML-SNAP)

    foreach(PKG ${ALL_PACKAGES})
      set(PKG_${PKG} ON CACHE BOOL "" FORCE)
    endforeach()

    pwd
cat <<EOT > ../cmake/presets/basic_milady.cmake
    # preset that turns on just a few, frequently used packages
    # this will be compiled quickly and handle a lot of common inputs.

    set(ALL_PACKAGES KSPACE MANYBODY MOLECULE RIGID EXTRA-FIX REPLICA ML-MILADY ML-SNAP)

    foreach(PKG ${ALL_PACKAGES})
      set(PKG_${PKG} ON CACHE BOOL "" FORCE)
    endforeach()
EOT
   
cat $MLD_LAMMPS_ROODIR/$TRG7/cmake/presets/basic_milady.cmake


cmake -C ../lammps_options.cmake -D PKG_MISC=on -D PKG_ML-MILADY=on -D PKG_ML-SNAP=on \
      -D CMAKE_INSTALL_PREFIX=$MLD_LAMMPS_ROODIR/$TRG7/install \
      -D CMAKE_C_COMPILER=$CCOMP \
      -D CMAKE_CXX_COMPILER=$CXXCOMP \
      -D CMAKE_Fortran_COMPILER=$FORTRANCOMP \
      -D ENABLE_TESTING=on \
      ../cmake
   
  
# cmake -C ../lammps_options.cmake -D PKG_MISC=on -D PKG_ML-MILADY=on -D PKG_ML-SNAP=on ../cmake
# cmake -C ../cmake/presets/basic_milady.cmake -D PKG_MISC=on -D PKG_ML-MILADY=on -D PKG_ML-SNAP=on ../cmake

cmake-gui ../cmake

cmake --build . --target all
cmake --build . --target install
cmake --build . --target test

# OK
cmake --build . --target test
    Running tests...
    
        Start 485: ImproperStyle:umbrella
485/486 Test #485: ImproperStyle:umbrella .............................   Passed    0.89 sec
        Start 486: ImproperStyle:zero
486/486 Test #486: ImproperStyle:zero .................................   Passed    0.91 sec

99% tests passed, 5 tests failed out of 486

Label Time Summary:
generated    =   0.88 sec*proc (1 test)
slow         =  58.17 sec*proc (63 tests)
unstable     =  29.44 sec*proc (32 tests)

Total Test time (real) = 424.31 sec

The following tests FAILED:
	227 - AtomicPairStyle:edip (Failed)
	235 - AtomicPairStyle:meam_spline (Failed)
	236 - AtomicPairStyle:meam_sw_spline (Failed)
	260 - ManybodyPairStyle:edip_multi (Failed)
	269 - ManybodyPairStyle:lcbop (Failed)
Errors while running CTest
Output from these tests are in: /home/catA/wambeke/MLD_LAMMPS/milady_lammps_cmake/build/Testing/Temporary/LastTest.log
Use "--rerun-failed --output-on-failure" to r
```


### Compilation ml_cv intel

on is246206

```
export MLD_LAMMPS_ROODIR=~/MLD_LAMMPS
export TRG6=ml_cv1/lammps
cd $MLD_LAMMPS_ROODIR/$TRG6


rm -rf build ; mkdir build ; cd build
bash

# on is246206 compile milady and run ctest and run python tests
# rename and modification and use at your own risk

# bash
# export PATH=/export/home/catA/wambeke/MLD/MILADY/scripts:${PATH}
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_intel_is246206.bash
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_intel_is246206.bash | grep -e '-- ..' | grep MLD_MKL_LIB
# f_stuff_intel_is246206

cmake --version
  cmake version 3.20.5


export LOC_ROODIR="/export/home/catA/wambeke"
  export MLD_OAP_ROODIR="/export/home/catA/intel"
  export MLD_OAP_INSDIR="${MLD_OAP_ROODIR}/oneapi"
  
export CONFIG_SETVARS=../config_setvars_oneapi.tmp
cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
mpi=latest
itac=latest
EOT

export CONFIG_SETVARS=../config_setvars_oneapi.tmp
cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
itac=latest
EOT

cat ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
itac=latest

source ${MLD_OAP_INSDIR}/setvars.sh --config=${CONFIG_SETVARS}
 
:: initializing oneAPI environment ...
   bash: BASH_VERSION = 5.1.0(1)-release
   args: Using "$@" for setvars.sh arguments: --config=../config_setvars_oneapi.tmp
:: compiler -- latest
:: mkl -- latest
:: mpi -- latest
:: itac -- latest
:: oneAPI environment initialized ::

CMAKE_C_COMPILER C compiler to be used for compilation (default: system specific, gcc on Linux)
CMAKE_CXX_COMPILER C++ compiler to be used for compilation (default: system specific, g++ on Linux)
CMAKE_Fortran_COMPILER Fortran compiler to be used for compilation (default: system specific, gfortran on Linux)

export CCOMP=$(which icc)
export CXXCOMP=$(which icpc)
export FORTRANCOMP=$(which ifort)
  #export FC=$(which ifort)
  #export CC=$(which icc)
  
export MKL_ROOT=${MKLROOT}

icc -help | head
icpc -help | head
ifort -help | head

export TRG6=ml_cv1/lammps

cd $MLD_LAMMPS_ROODIR/$TRG6
rm -rf build ; mkdir build ; cd build

cat ../cmake/presets/basic.cmake
    # preset that turns on just a few, frequently used packages
    # this will be compiled quickly and handle a lot of common inputs.

    set(ALL_PACKAGES KSPACE MANYBODY MOLECULE RIGID)

    foreach(PKG ${ALL_PACKAGES})
      set(PKG_${PKG} ON CACHE BOOL "" FORCE)
    endforeach()
    
cat <<EOT > ../cmake/presets/basic_milady.cmake    
    # preset that turns on just a few, frequently used packages
    # this will be compiled quickly and handle a lot of common inputs.

    set(ALL_PACKAGES EXTRA-FIX KSPACE MANYBODY MISC ML-MILADY ML-SNAP MOLECULE REPLICA RIGID USER_MILADY)

    foreach(PKG ${ALL_PACKAGES})
      message(STATUS "!!!! preset on - ${PKG}")
      set(PKG_${PKG} ON CACHE BOOL "on" FORCE)
    endforeach()
EOT
   
cat $MLD_LAMMPS_ROODIR/$TRG6/cmake/presets/basic_milady.cmake
```

#### compile intel or gnu with mpi or openmpi

```
cd /home/catA/wambeke/MLD_LAMMPS/ml_cv1/lammps
rm -rf build ; mkdir build ; cd build

# new all in one
 
module avail
  ---------- /usr/share/Modules/modulefiles 
  dot  module-git  module-info  modules  null  use.own  
  ---------- /usr/share/modulefiles 
  mp-x86_64  mpi/mpich-x86_64  mpi/openmpi-x86_64 

# for gnu is useful to get mpi, else is OK also and 293 tests are moooore speedy (20s vs 1035s!)
bash 
module load mpi/openmpi-x86_64

# new all in one
cmake -C ../cmake/presets/basic_milady.cmake \
  -D CMAKE_VERBOSE_MAKEFILE=off \
  -D CMAKE_INSTALL_PREFIX=$MLD_LAMMPS_ROODIR/$TRG6/install \
  -D ENABLE_TESTING=on \
  ../cmake
  
# intel  
# for intel is useful to get mpi, 296 tests are moooore speedy (452s vs 1035s!) 
cmake -C ../cmake/presets/basic_milady.cmake \
  -D CMAKE_VERBOSE_MAKEFILE=off \
  -D CMAKE_INSTALL_PREFIX=$MLD_LAMMPS_ROODIR/$TRG6/install \
  -D CMAKE_C_COMPILER=$CCOMP \
  -D CMAKE_CXX_COMPILER=$CXXCOMP \
  -D CMAKE_Fortran_COMPILER=$FORTRANCOMP \
  -D ENABLE_TESTING=on \
  ../cmake
   
cmake-gui ../cmake

# --parallel 4
# nproc --> 20
# mate-system-monitor &

cmake --build . --target all --parallel 15
cmake --build . --target install
cmake --build . --target test

# OK
cmake --build . --target test
    Running tests...
    Test project /home/catA/wambeke/MLD_LAMMPS/ml_cv1/lammps/build
            Start   1: RunLammps
      1/495 Test   #1: RunLammps ..........................................   Passed    0.88 sec

            Start 490: ImproperStyle:hybrid
    490/495 Test #490: ImproperStyle:hybrid ...............................   Passed    0.88 sec
            Start 491: ImproperStyle:inversion_harm
            nic
    491/495 Test #491: ImproperStyle:inversion_harmonic ...................   Passed    0.87 sec
            Start 492: ImproperStyle:ring
    492/495 Test #492: ImproperStyle:ring .................................   Passed    0.87 sec
            Start 493: ImproperStyle:sqdistharm
    493/495 Test #493: ImproperStyle:sqdistharm ...........................   Passed    0.88 sec
            Start 494: ImproperStyle:umbrella
    494/495 Test #494: ImproperStyle:umbrella .............................   Passed    0.87 sec
            Start 495: ImproperStyle:zero
    495/495 Test #495: ImproperStyle:zero .................................   Passed    0.88 sec

    100% tests passed, 0 tests failed out of 495

    Label Time Summary:
    noWindows    =   1.78 sec*proc (2 tests)
    slow         =  55.88 sec*proc (64 tests)
    unstable     =  29.67 sec*proc (34 tests)

    Total Test time (real) = 426.61 sec
    
The following tests FAILED:
	 32 - ComputeGlobal (Failed)
	222 - AtomicPairStyle:edip (Failed)
	230 - AtomicPairStyle:meam_spline (Failed)
	231 - AtomicPairStyle:meam_sw_spline (Failed)
	261 - ManybodyPairStyle:edip_multi (Failed)
	270 - ManybodyPairStyle:lcbop (Failed)
	350 - KSpaceStyle:ewald_tri (Failed)
	388 - FixTimestep:addtorque_const (Failed)
	407 - FixTimestep:nph (Failed)
	410 - FixTimestep:npt_iso (Failed)
	436 - FixTimestep:rigid_npt_small (Failed)
	458 - FixTimestep:temp_csld (Failed)
	481 - DihedralStyle:table_linear (Failed)
	482 - DihedralStyle:table_spline (Failed)
	
ctest --rerun-failed --output-on-failure
Test project /home/catA/wambeke/MLD_LAMMPS/ml_cv1/lammps/build
      Start  32: ComputeGlobal
 1/14 Test  #32: ComputeGlobal ....................***Failed    0.90 sec
 etc...

```


### Compilation lammps_stable_cmake

in lammps/doc/Manual.pdf

```
8.6 Tutorials howto
8.6.1 Using CMake with LAMMPS tutorial

Some common CMake variables
CMAKE_INSTALL_PREFIX root directory of install location for make install (default: $HOME/.local)
CMAKE_BUILD_TYPE controls compilation options: one of RelWithDebInfo (default), Release, Debug, MinSizeRel
CMAKE_VERBOSE_MAKEFILE if set to on echo commands while executing during build (default: off)


Example: cmake --build . --target all or make all. 
The following abstract targets are available:
all, lammps, doc, install, test, clean

build “everything” (default)
build the LAMMPS library and executable
build the HTML documentation (if configured)
install all target files into folders in CMAKE_INSTALL_PREFIX
run some tests (if configured with -D ENABLE_TESTING=on)
remove all generated files
```

Compilation

```
export TRG5=lammps_stable_compile_cmake
export TRG5=lammps_stable_compile_cmake_intel

export MAKEFLAGS=-j5
unset MAKEFLAGS

cd $MLD_LAMMPS_ROODIR
rm -rf $TRG5
cp -r lammps_stable $TRG5
cd $TRG5
git his1
    88c8b6ec6f 2022-11-03 | Merge pull request #3460 from akohlmey/maintenance-2022-06-23 
    (HEAD -> stable, tag: stable_23Jun2022_update2, tag: patch_23Jun2022_update2, origin/stable) <Axel Kohlmeyer>
rm -rf build
mkdir build
cd build

cmake ../cmake -D CMAKE_INSTALL_PREFIX=$MLD_LAMMPS_ROODIR/$TRG5/install -D ENABLE_TESTING=on
cmake -C ../cmake/presets/basic.cmake -D PKG_MISC=on ../cmake
# cmake --build --parallel 4 . --target all  # New from cmake v3.12 # --> parallel KO
cmake --build . --target all
cmake --build . --target install
cmake --build . --target test

  --> Ok sauf pour test_pair_style référence à « _intel_fast_memcpy » non définie
  --> TODO compiler with intel...
  
cmake --build . --target test

    Running tests...
    Test project /home/catA/wambeke/MLD_LAMMPS/lammps_stable_compile_cmake/build
            Start   1: RunLammps
      1/496 Test   #1: RunLammps ..........................................   Passed    1.74 sec
            Start   2: HelpMessage
      2/496 Test   #2: HelpMessage ........................................   Passed    1.49 sec
      
...

            Start  45: FortranOpen
     45/496 Test  #45: FortranOpen ........................................   Passed    1.62 sec
            Start  46: FortranCommands
     46/496 Test  #46: FortranCommands ....................................   Passed    1.62 sec
            Start  47: ErrorStats
     47/496 Test  #47: ErrorStats .........................................   Passed    0.00 sec
            Start  48: MolPairStyle:beck
    Could not find executable /home/catA/wambeke/MLD_LAMMPS/lammps_stable_compile_cmake/build/test_pair_style
     9% tests passed, 449 tests failed out of 496
     

 
```

### Compilation lammps_stable_cmake_intel


```
@is246206.../lammps_stable_compile_cmake_intel(stable)>

rm -rf build
mkdir build
cd build
@is246206.../build(stable)>bash
@is246206.../build(stable)> 
# on is246206 compile milady and run ctest and run python tests
# rename and modification and use at your own risk

# bash
# export PATH=/export/home/catA/wambeke/MLD/MILADY/scripts:${PATH}
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_intel_is246206.bash
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_intel_is246206.bash | grep -e '-- ..' | grep MLD_MKL_LIB
# f_stuff_intel_is246206

function in_red () {
  # Red bold
  echo -e "\e[31m\e[1m"$@"\e[0m"
  
function in_red_pause () {
  # Red bold
  echo -e "\e[31m\e[1m"$@"\e[0m"
  f_pause
}

function f_pause {
  read -p 'Continue? (ctrl-c to stop)' rep
}

function in_green () {
  # green bold
  echo -e "\e[32m\e[1m"$@"\e[0m"
}
@is246206.../build(stable)> cmake --version
cmake version 3.20.5

CMake suite maintained and supported by Kitware (kitware.com/cmake).
@is246206.../build(stable)> export LOC_ROODIR="/export/home/catA/wambeke"
  export MLD_OAP_ROODIR="/export/home/catA/intel"
  export MLD_OAP_INSDIR="${MLD_OAP_ROODIR}/oneapi"
@is246206.../build(stable)> export CONFIG_SETVARS=../config_setvars_oneapi.tmp
@is246206.../build(stable)> cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
mpi=latest
itac=latest
EOT
@is246206.../build(stable)> cat ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
mpi=latest
itac=latest
@is246206.../build(stable)> source ${MLD_OAP_INSDIR}/setvars.sh --config=${CONFIG_SETVARS}
 
:: initializing oneAPI environment ...
   bash: BASH_VERSION = 5.1.0(1)-release
   args: Using "$@" for setvars.sh arguments: --config=../config_setvars_oneapi.tmp
:: compiler -- latest
:: mkl -- latest
:: mpi -- latest
:: itac -- latest
:: oneAPI environment initialized ::

CMAKE_C_COMPILER C compiler to be used for compilation (default: system specific, gcc on Linux)
CMAKE_CXX_COMPILER C++ compiler to be used for compilation (default: system specific, g++ on Linux)
CMAKE_Fortran_COMPILER Fortran compiler to be used for compilation (default: system specific, gfortran on Linux)

export CCOMP=$(which icc)
export CXXCOMP=$(which icpc)
export FORTRANCOMP=$(which ifort)
  #export FC=$(which ifort)
  #export CC=$(which icc)
  
export MKL_ROOT=${MKLROOT}

icc -help | head
icpc -help | head
ifort -help | head

cmake ../cmake -D CMAKE_INSTALL_PREFIX=$MLD_LAMMPS_ROODIR/$TRG5/install \
               -D CMAKE_C_COMPILER=$CCOMP \
               -D CMAKE_CXX_COMPILER=$CXXCOMP \
               -D CMAKE_Fortran_COMPILER=$FORTRANCOMP \
               -D ENABLE_TESTING=on
cmake -C ../cmake/presets/basic.cmake -D PKG_MISC=on ../cmake

# avoid  'export MAKEFLAGS=-j8' or 'cmake --build . --parallel 4 --target all' # --> parallel KO

cmake --build . --target all
cmake --build . --target install
cmake --build . --target test

# presque OK

make --build . --target test
  Running tests...
  Test project /home/catA/wambeke/MLD_LAMMPS/lammps_stable_compile_cmake_intel/build
          Start   1: RunLammps
    1/496 Test   #1: RunLammps ..........................................   Passed    0.90 sec
          Start   2: HelpMessage
    2/496 Test   #2: HelpMessage ........................................   Passed    0.90 sec
          Start   3: InvalidFlag

  494/496 Test #494: ImproperStyle:sqdistharm ...........................   Passed    0.88 sec
          Start 495: ImproperStyle:umbrella
  495/496 Test #495: ImproperStyle:umbrella .............................   Passed    0.88 sec
          Start 496: ImproperStyle:zero
  496/496 Test #496: ImproperStyle:zero .................................   Passed    0.87 sec

  98% tests passed, 12 tests failed out of 496

  Label Time Summary:
  noWindows    =   1.74 sec*proc (2 tests)
  slow         =  67.37 sec*proc (64 tests)
  unstable     =  31.27 sec*proc (34 tests)

  Total Test time (real) = 445.16 sec

  The following tests FAILED:
	   32 - ComputeGlobal (Failed)
	  222 - AtomicPairStyle:edip (Failed)
	  230 - AtomicPairStyle:meam_spline (Failed)
	  231 - AtomicPairStyle:meam_sw_spline (Failed)
	  261 - ManybodyPairStyle:edip_multi (Failed)
	  270 - ManybodyPairStyle:lcbop (Failed)
	  350 - KSpaceStyle:ewald_tri (Failed)
	  407 - FixTimestep:nph (Failed)
	  410 - FixTimestep:npt_iso (Failed)
	  436 - FixTimestep:rigid_npt_small (Failed)
	  481 - DihedralStyle:table_linear (Failed)
	  482 - DihedralStyle:table_spline (Failed)
  Errors while running CTest


```


### compilation milady/lammps


```
cat cmake/presets/basic.cmake 
    # preset that turns on just a few, frequently used packages
    # this will be compiled quickly and handle a lot of common inputs.

    set(ALL_PACKAGES KSPACE MANYBODY MOLECULE RIGID)

    foreach(PKG ${ALL_PACKAGES})
      set(PKG_${PKG} ON CACHE BOOL "" FORCE)
    endforeach()
    
    
cat milady_lammps/lammps_options.cmake
    # FOR LAMMPS CMAKE BUILD

    # set installation location
    set(CMAKE_INSTALL_PREFIX "$ENV{PREFIX}")

    # enforce c++11 standards
    set(CCFLAGS -g -O3 -std=c++11)

    # compile a binary and shared library
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

    # allow error messages (very useful)
    set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)

    # minimal packages to run example (MANYBODY, EXTRA-FIX) and
    # generate new pathway (REPLICA for "fix neb")
    # Also include some "ML" potentials (PACE and SNAP)
    set(ALL_PACKAGES MANYBODY EXTRA-FIX REPLICA ML-MILADY ML-SNAP)

    foreach(PKG ${ALL_PACKAGES})
      set(PKG_${PKG} ON CACHE BOOL "" FORCE)
    endforeach()
    
    

```
