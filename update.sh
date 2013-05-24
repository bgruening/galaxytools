cd src
jekyll
cd ..
cp src/_site/* ./ -r
git add .
git commit -m 'website updates'
git push
