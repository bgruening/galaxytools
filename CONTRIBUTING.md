# Contributing

This document describes how to contribute to this repository. Pull
requests containing bug fixes, updates, and extensions to the existing
tools and tool suites in this repository will be considered for
inclusion.

## How to Contribute

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Make sure you have git [installed](https://help.github.com/articles/set-up-git)
* Fork the repository on [GitHub](https://github.com/bgruening/galaxytools/fork)
* Make the desired modifications - consider using a [feature branch](https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches).
* Try to stick to the [IUC standards](http://galaxy-iuc-standards.readthedocs.org/en/latest/) whenever you can
* Make sure you have added the necessary tests for your changes and they pass.
* Open a [pull request](https://help.github.com/articles/using-pull-requests)
  with these changes.

## Setting up GitHub Environment for Deployment (Maintainers Only)

If you are a regular contributor working on a fork, you do **not** need to set this up. The deployment
jobs will only run on the main upstream repository when code is merged to the `master` or `main` branch. 
All other CI checks (linting, testing) will work normally on your fork without any special configuration.

### Background

This repository uses GitHub Environments instead of repository-level secrets
for deployment. This provides better security through deployment protection rules and clearer
separation of deployment credentials from general repository secrets.

**This environment is only configured on the main upstream repository and should not be replicated on forks.**

### Setup Instructions

1. **Create the Environment**
   - Go to your repository: `Settings` > `Environments`
   - Click `New environment`
   - Name it: `toolshed-deployment`

2. **Add Required Secrets to the Environment**
   
   **IMPORTANT**: All secrets must be added to the `toolshed-deployment` **environment**,
   not as repository-level secrets (Settings > Secrets and variables > Actions).
   
   Go to: `Settings` > `Environments` > `toolshed-deployment` > `Add Secret`
   
   Add each of the following secrets:
   
   - **TTS_API_KEY**: Your API key for TestToolShed
   - **TS_API_KEY**: Your API key for ToolShed
   - **PAT**: Personal Access Token with `repo` scope
