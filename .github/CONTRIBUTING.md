# Contributing

Thank you for your interest in contributing to WHAM! 

We welcome simple documentation and bug fixes as well as more substantive code contributions, such as new features. Our current ideas for improving WHAM are on the [to do list](https://github.com/timjmiller/wham/wiki/To-do-list). If you would like to tackle something larger than a bug fix, please file an [Issue](https://github.com/timjmiller/wham/issues) to discuss and make sure we are not working on it already.

We reserve the **master** branch for *stable, tagged [releases](https://github.com/timjmiller/wham/releases)*. The **devel** branch contains *unstable but tested* updates. We create additional, temporary branches to develop and test larger updates before merging them into devel. While working on a temporary branch, we semi-regularly pull updates from devel to deal with merge conflicts as soon as possible.

If you wish to contribute, please use the NOAA Fisheries Toolbox suggested workflow, summarized below. [This page](https://github.com/nmfs-fish-tools/Resources/blob/master/CONTRIBUTING.md) has further guidance and resources. The "fork-and-branch" workflow we follow is described in more detail [here](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/).

## Workflow

1. File an [Issue](https://github.com/timjmiller/wham/issues) describing the fix/feature you are planning.
2. Fork wham (create your own copy on GitHub by clicking "Fork", top right of wham page)
3. Clone your fork (download your copy of wham so you can modify it locally, click green "Code" button)
4. Add the original wham repo as *upstream* so you can get future changes from devel (see [here](https://docs.github.com/en/github/collaborating-with-pull-requests/working-with-forks/configuring-a-remote-for-a-fork)).
    
    ```
    git remote add upstream https://github.com/timjmiller/wham.git
    ``` 
5. Checkout devel.
    
    ```
    git checkout devel
    ```
6. Confirm you can fetch and merge changes from the original wham repo (see [here](https://docs.github.com/en/github/collaborating-with-pull-requests/working-with-forks/syncing-a-fork)).
    
    ```
    git fetch upstream
    git merge upstream/devel
    ``` 
7. Create a new branch for your change with a meaningful name, e.g. `selftest` to add a function to simulate data from a fit model and refit.
    
    ```
    git checkout -b <myfeature>
    ```
8. Confirm you can run and pass the automated tests
    
    ```
    # install.packages("here")
    library(here)
    wham.source.dir <- here()
    
    # confirm that the .cpp compiles and wham can be loaded
    devtools::load_all(path = wham.source.dir)
    
    # run all tests
    devtools::test(pkg = wham.source.dir)
    ```
9. Work on your feature branch, making changes to your local wham files. Periodically merge in updates from the upstream devel branch (#6 above). Try to maintain the current wham coding style.
10. Add documentation, test(s), and a vignette if applicable.
11. When your feature is finished, rerun the tests and verify that they pass (#8). Helpful scripts for testing and debugging are in the `inst/contribute/` folder. 
12. Confirm your branch is up-to-date with the upstream devel.
13. Commit your changes, and push your local branch to your forked repo (origin, not upstream)
    
    ```
    git commit -m "informative message"
    git push origin <myfeature>
    ```
14. Submit your pull request, following these [best practices](https://www.atlassian.com/blog/git/written-unwritten-guide-pull-requests).
15. Clean up / delete your feature branch and fork.
