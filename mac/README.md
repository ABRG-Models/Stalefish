To create a Mac release (these instructions are mainly for me [Seb]):

Recompile **and install** Stalefish after git-pulling Stalefish and, if necessary, morphologica.

```
cd models/Stalefish
git pull
cd build
make
sudo make install
cd ..
```

In the mac directory, run bundle.sh to make the Mac application bundle. This picks out the binaries that you just built and makes the bundle directory tree. It depends on the binaries having been installed exactly like I did.

```
cd mac
./bundle.sh
```

Make a .dmg file. That means following the recipe from https://medium.com/@mattholt/packaging-a-go-application-for-macos-f7084b00f6b5

Basically, you place the Stalefish.app directory and an Applications link into a .dmg disk image.
