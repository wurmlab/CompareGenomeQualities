FROM continuumio/miniconda3:4.9.2

# Install curl and latex through apt.
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl texlive-latex-base texlive-latex-recommended \
    texlive-fonts-recommended lmodern cm-super \
    && rm -rf /var/lib/apt/lists/*

# Copy repo to image.
COPY . /cmg

# Install the tools we need.
RUN /cmg/conda_install.sh

# Create empty ~/.parallel/ignored_vars and ~/.parallel/will-cite files so
# parallel will pass all env variables and not show the citation message.
# And yes, I agree to cite parallel.
RUN mkdir -p ~/.parallel \
&& touch ~/.parallel/will-cite \
&& touch ~/.parallel/ignored_vars

ENTRYPOINT ["/cmg/compare-genome-qualities.sh"]
WORKDIR /mnt
