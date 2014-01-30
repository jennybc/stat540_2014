Git and GitHub
==============

> Git allows groups of people to work on the same documents (often code)
at the same time, and without stepping on each other's toes. It's a
distributed version control system.

(cribbed from [tryGit][trygit])

Git
===

Install the git command line program
------------------------------------

### Mac

+ [Install Xcode](macappstores://itunes.apple.com/us/app/xcode/id497799835)
  from the App Store
+ Open Xcode, and select
  *Xcode -> Preferences -> Downloads -> Command Line Tools -> Install*
+ Run `git --version` at the terminal prompt to verify that git's installed

### Windows

Install [Git for Windows][msysgit] (also known as msysgit) to install
the git command line program needed by RStudio.

Alternatively, install [GitHub for Windows][githubwindows], which
includes both a GitHub GUI and the git command line program.

[msysgit]: http://msysgit.github.io/
[githubwindows]: http://windows.github.com/

### Ubuntu or Debian Linux

```sh
sudo apt-get install git
```

### Fedora or RedHat Linux

```sh
sudo yum install git
```

Configure the git command line program
--------------------------------------

Follow the [instructions on GitHub][ghsetup] to configure git.

[ghsetup]: https://help.github.com/articles/set-up-git

Learning to use git at the command line
---------------------------------------

We won't be exploring git at the command line in this course, but it
is a rather useful skill, particularly to get yourself out of sticky
situations involving conflicts and merges.

+ [tryGit][trygit] is a fantastic interactive tutorial for learning to
  use git at the command line.

[trygit]: http://try.github.io/

Using git with RStudio
----------------------

Read the RStudio documentation on using git with RStudio:

+ [Using Version Control with RStudio][rstudiogit]

+ [Create a repository][ghcreate] on GitHub, and select the
  checkbox labeled *Initialize this repository with a README*
+ Open the repository in your web browser
+ [Copy the repository URL][ghclone] of your project to the clipboard by clicking the *Copy
  to clipboard* icon next to *HTTPS clone URL* on the right-hand side
  of the project page
+ Follow the [RStudio documentation][rstudiogit] for
  *Creating a new project based on a Git or Subversion repository*
+ It's also possible to
  + create a new RStudio project from scratch and request that RStudio make it a git repo upon creation
  + or add version control to an existing RStudio project
  + and then import that git repository into GitHub

[rstudiogit]: http://www.rstudio.com/ide/docs/version_control/overview
[ghcreate]: https://help.github.com/articles/create-a-repo
[ghclone]: https://help.github.com/articles/fetching-a-remote

GitHub password
---------------

You may see this error message from RStudio:

```
error: unable to read askpass response from 'rpostback-askpass'
fatal: could not read Username for 'https://github.com': Device not configured
```

There's a few ways to go about fixing the password issue.

### Password Caching

The easiest is to tell git to cache your password for the next four
months (about ten million seconds) by opening a terminal and running
the following command:

```sh
git config --global credential.helper 'cache --timeout=10000000'
```

You'll then have to force git to ask for your password once at the
command line. You can do this by cloning your GitHub repository and
then pushing it back to GitHub. In the following commands, substitute
USERNAME for your GitHub user name, and REPOSITORY for the name of
your stat540 repository.

```sh
git clone https://github.com/USERNAME/REPOSITORY
cd REPOSITORY
git push
```

These instructions are on the GitHub help page under
[Password Caching][githubpassword].

### SSH RSA Key

Another option is to set up git, GitHub and RStudio to use an SSH RSA key.
In RStudio, click *Preferences -> Git/SVN -> SSH RSA Key -> Create RSA
Key* then *View public key* and then follow these instructions to add
[Add your SSH key to GitHub][githubsshkey]. You may need to restart
RStudio before the configuration takes effect.

[githubpassword]: https://help.github.com/articles/set-up-git#platform-linux
[githubsshkey]: https://help.github.com/articles/generating-ssh-keys#step-3-add-your-ssh-key-to-github

### Mac OS Keychain Credential Helper

On Mac OS, you have one more option. You can configure git to use the
Mac OS keychain helper `git-credential-osxkeychain` to store your
password. Folllow the GitHub Help instructions for
[Password caching][githubosxkeychain].

[githubosxkeychain]: https://help.github.com/articles/set-up-git#platform-mac

Install a GUI for git
---------------------

We're recommending that you use RStudio to manage your files in
GitHub, which keeps everything all in one nice tidy package. That
being said, there are other pretty GUIs out there for git and GitHub.

Atlassian SourceTree is a nice complement to RStudio's somewhat
limited Git GUI. Some useful features of SourceTree are:

+ shows visual evidence when your local repo is not current
  and there are remote commits that you need to pull
+ launch external merge tools to resolve conflicts
+ configure [remote repositories][gitscmremotes] and
  [remote tracking branches][gitscmremotebranches]

[gitscmremotes]: http://git-scm.com/book/en/Git-Basics-Working-with-Remotes
[gitscmremotebranches]: http://git-scm.com/book/en/Git-Branching-Remote-Branches

### Mac

+ [Atlassian SourceTree](http://www.sourcetreeapp.com/)
+ [GitHub for Mac](http://mac.github.com/)
+ [GitX](http://gitx.frim.nl/)

### Windows

+ [Atlassian SourceTree](http://www.sourcetreeapp.com/)
+ [GitHub for Windows](http://windows.github.com/)

### Linux

+ [gitg](http://live.gnome.org/Gitg)
+ gitk comes built-in with git

### More GUIs

There's more git GUIs listed on
[git-scm.com](http://git-scm.com/downloads/guis).

GitHub
======

Getting Started with GitHub
---------------------------

+ Register an [educational GitHub account](http://GitHub.com/edu)
+ Upload your photo to [Gravatar](http://gravatar.com), if you're comfortable with that
+ Create a public GitHub repository named `stat540-2014-lastname-firstname` for lab work
+ Create a private GitHub repository named `stat540-2014-lastname-firstname-hw` for assignments
+ Add the profs and TAs as collaborators to your private GitHub repository
  - [Jenny Bryan](https://github.com/jennybc) @jennybc on GitHub
  - [Shaun Jackman](https://github.com/sjackman) @sjackman on GitHub
  - [Gloria Li](https://github.com/gloriali) @gloriali on GitHub

GitHub Documentation
--------------------

GitHub has great documentation on both git and GitHub.

+ [GitHub Help](https://help.github.com/)
+ [Set Up Git][ghsetup]
+ [Create A Repository][ghcreate]
+ [Fetching a remote][ghclone]

Further Reading
===============

+ [An introduction to Git/Github](http://kbroman.github.io/github_tutorial/)
  by Karl Broman, aimed at stats / data science types
+ Ram, K 2013. Git can facilitate greater reproducibility and
  increased transparency in science. Source Code for Biology and
  Medicine 2013 8:7. Go to the
  [associated github repo](https://github.com/karthikram/smb_git)
  to get the PDF (link at bottom of README) and to see a great example
  of how someone managed the process of writing a manuscript with
  git(hub).
+ Blog post [Version control for scientific research](http://blogs.biomedcentral.com/bmcblog/2013/02/28/version-control-for-scientific-research/)
  by Karthik Ram and Titus Brown on the BioMed Central blog February 28
+ Blog post [Getting Started With a GitHub Repository](http://chronicle.com/blogs/profhacker/getting-started-with-a-github-repository)
  from ProfHacker on chronicle.com looks helpful
+ Blog post from The Molecular Ecologist on
  [using GitHub with R and RStudio](http://www.molecularecologist.com/2013/11/using-github-with-r-and-rstudio/)

## [PHD Comics - "FINAL".doc](http://www.phdcomics.com/comics/archive.php?comicid=1531)

!["FINAL".doc](http://www.phdcomics.com/comics/archive/phd101212s.gif)

## [XKCD - Git Commit](http://xkcd.com/1296/)

![Git Commit](http://imgs.xkcd.com/comics/git_commit.png)
