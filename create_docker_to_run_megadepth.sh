#to be run after ./run_hbb_build_container.sh (produces actual megadepth binary referenced here)
mkdir -p docker.run.build
rsync -av megadepth docker.run.build/
rsync -av Dockerfile.run docker.run.build/Dockerfile
pushd docker.run.build/
VER=$(cat ../VERSION)
docker build --tag quay.io/broadsword/megadepth:${VER} --tag quay.io/broadsword/megadepth:latest .
popd

docker push quay.io/broadsword/megadepth:${VER}
docker push quay.io/broadsword/megadepth:latest
