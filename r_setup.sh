#!/usr/bin/env bash
set -e

echo "🔄 Updating Ubuntu..."
sudo apt update && sudo apt upgrade -y

echo "📦 Installing system dependencies for R packages..."
sudo apt install -y \
build-essential cmake gfortran r-base-dev cargo \
libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
libblas-dev liblapack-dev libatlas-base-dev \
libpng-dev libjpeg-dev libtiff5-dev libcairo2-dev libfreetype6-dev \
libfontconfig1-dev libharfbuzz-dev libfribidi-dev libxt-dev \
libudunits2-dev libgsl-dev \
libmysqlclient-dev libpq-dev \
librsvg2-dev libgeos-dev libproj-dev libgdal-dev \
ffmpeg libpoppler-cpp-dev tesseract-ocr libtesseract-dev \
imagemagick libmagick++-dev

echo "✅ System libraries installed."

echo "📦 Installing pak (for fast parallel R installs)..."
Rscript -e 'install.packages("pak", repos="https://cloud.r-project.org")'

echo "📦 Installing common R packages..."
Rscript -e 'pak::pak(c(
  "tidyverse", "data.table", "MatrixModels", "quantreg", "car", "alr4",
  "rio", "caret", "glmnet", "lme4", "nlme", "survival",
  "magick", "av", "gifski", "pdftools", "tesseract",
  "sf", "units", "gsl", "RMySQL", "RPostgreSQL", "rsvg",
  "plotly", "RColorBrewer", "ggthemes", "devtools", "remotes"
))'

echo "🎉 R environment setup complete!"


