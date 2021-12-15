# ATL1415
Code used to generate NASA's ATL14 and ATL15 products for the ICESat-2 project

This repository contains code originally published under the smithB/surfaceChange repository.  The switch to the new repository is intended to mark release 1.0 of the products.  Work on this project has involved major contributions from git users bpjelley, tsutterley, and suzanne64.

##ATL14/15 products

ATL14 and ATL15 provide surface height and surface-height data for ice sheets in the Antarctic and the Arctic, as measured by NASA's ICESat'2 project.  ATL14 is a reference digital elevation model (DEM) with a time stamp of 2020.0.  ATL15 provides surface-height change relative to ATL14.  

The products are described in an Algorithm Theoretical Basis Document that is available at the National Snow and Ice Data center (see https://nsidc.org/data/ATL14 and https://nsidc.org/data/ATL15)

##Algorithm demonstrations

Jupyter notebooks that demonstrate the algorithm are in the /notebooks directory.  A good place to start is: [Insert]




### Contributing Code or Examples

We follow a standard Forking Workflow for code changes and additions.
Submitted code goes through the pull request process for evaluation and comments.

#### General Guidelines

- Make each pull request as small and simple as possible
- Commit messages should be clear and describe the changes
- Larger changes should be broken down into their basic components and integrated separately
- If possible, bug fixes should be their own pull requests with an associated [GitHub issue](https://github.com/SmithB/surfaceChange/issues)
- Write a descriptive pull request message with a clear title

#### Steps to Contribute

1. Fork the repository to your personal GitHub account by clicking the "Fork" button on the project [main page](https://github.com/SmithB/surfaceChange).  This creates your own server-side copy of the repository.
2. Create a work environment to make your changes by cloning to your local system.
3. Add your fork as the `origin` remote and the original project repository as the `upstream` remote.  While this step isn't a necessary, it allows you to keep your fork up to date in the future.
4. Create a new branch to do your work.
5. Make your changes on the new branch and add yourself to the list of [project contributors](./CONTRIBUTORS.md).
6. Push your work to GitHub under your fork of the project.
7. Submit a [Pull Request](https://github.com/SmithB/surfaceChange/pulls) from your forked branch to the project repository.
