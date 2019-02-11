.. _yade-github-label:

##############
Yade on GitLab
##############

************************************************
Fast checkout without GitLab account (read-only)
************************************************
 
Getting the source code without registering on GitLab can be done via a single command. It will not allow interactions with the remote repository, which you access the read-only way::

 git clone https://gitlab.com/yade-dev/trunk.git

************************
Using branches on GitLab
************************

Most usefull commands are below. For more details, see for instance http://gitref.org/index.html and https://help.github.com/articles/set-up-git

Setup
=====

1. Register on gitlab.com

2. Add your `SSH key <https://gitlab.com/profile/keys>`_ to GitLab

3. Set your username and email through terminal

   ::
 
      git config --global user.name "Firstname Lastname"
      git config --global user.email "your_email@youremail.com"
 
   You can check these settings with ``git config --list``.


4. `Fork <https://help.github.com/articles/fork-a-repo>`_ the repository:

   Click the “Fork” button on the `gitlab page <https://gitlab.com/yade-dev/trunk>`_
   
   .. note:: By default gitlab will try and compile the forked repository, and it will fail if you don't have runners attached to your account. To avoid receiving failure notifications go to repository settings (bottom of left panel->general->permissions) to turn of pipelines. 

5. Set Up Your Local Repo through terminal:

   ::

      git clone git@gitlab.com:username/trunk.git
   This creates a new folder, named trunk, that contains the whole code (make sure username is replaced by your GitLab name).

6. Configure remotes

   ::

    cd to/newly/created/folder
    git remote add upstream git@gitlab.com:yade-dev/trunk.git
    git remote update
    
   Now, your "trunk" folder is linked with two remote repositories both hosted on gitlab.com, the original trunk from yade-dev (called "upstream" after the last command) and the fork which resides in your personal account (called "origin" and always configured by default). Through appropriate commands explained below, you will be able to update your code to include changes commited by others, or to commit yourself changes that others can get.

Retrieving older Commits
========================

In case you want to work with, or compile, an older version of Yade which is not tagged, you can create your own (local) branch of the corresponding daily build. Look `here <https://answers.launchpad.net/yade/+question/235867>`_ for details.

Committing and updating 
========================

For those used to other version control systems, note that the commit mechanisms in Git significantly differs from that of `Bazaar <http://bazaar.canonical.com/en/>`_ or `SVN <https://subversion.apache.org/>`_. Therefore, don't expect to find a one-to-one command replacement. In some cases, however, the equivalent bazaar command is indicated below to ease the transition.

Inspecting changes
------------------

You may start by inspecting your changes with a few commands. For the "diff" command, it is convenient to copy from the output of "status" instead of typing the path to modified files. ::

   git status
   git diff path/to/modified/file.cpp

Pushing changes to remote repository
------------------------------------

1. Push to personnal repository

   After previous steps proceed to commit through terminal::

      git add path/to/new/file.cpp  #Version a newly created file
      git commit path/to/new_or_modified/file.cpp -m 'Commit message'  #stage (register) change in the local repository
      git push  #Push all commits to the remote branch
  
   The changes will be pushed to your personal fork.


2. Push to yade-dev

   Merging changes into yade-dev's master branch cannot be done directly with a push, only by merge request (see below). It is possible however to push changes to a new branch of yade-dev repository for members of that group. It is `currently <https://gitlab.com/gitlab-org/gitlab-ce/issues/23902>`_ the only way to have merge requests tested by the gitlab CI pipeline before being effectively merged. To push to a new yade-dev/branch::

      git push upstream yourBranch:newlyCreatedBranch #Push all commits to a new remote branch.

   in that command ``yourBranch`` is the local branch name on which you were working (possibly ``master``) and ``newlyCreatedBranch`` will be the name of that branch on the remote. Please choose a descriptive name as much as you can.

Requesting merge into yade-dev master branch
--------------------------------------------

If you have tested your changes and you are ready to merge them into yade-dev's master branch, you'll have to make a "merge request" (MR) from the gitlab.com interface (see the "+" button at the top of the repository webpage). Set source branch and target branch, from yade-dev/trunk/newlyCreatedBranch to yade-dev/trunk/master. The MR will trigger a `pipeline <https://gitlab.com/yade-dev/trunk/pipelines>`_ which includes compiling, running regression tests, and generating the documentation (the `newly built <https://yade-dev.gitlab.io/trunk>`_ documentation is accessible via settings->pages).
If the full pipeline succeeds the merge request can be merged into master branch.

.. note::
   In case of MR to yade-dev's master from another branch of yade-dev, the pipeline will use group runners attached to yade-dev (the group runners are kindly provided by `3SR <https://www.3sr-grenoble.fr/?lang=en>`_ and `UMS Gricad <https://gricad.univ-grenoble-alpes.fr/>`_).
   If the MR is from a branch of a forked repository (under personnal account) however, the pipeline needs runners available under the personnal account (check this with your local IT support). If you don't have access to gitlab runners pushing to a branch of yade-dev is mandatory (method 2 in previous section).

Alternatively, create a patch from your commit via::

 git format-patch origin  #create patch file in current folder)

and send to the developers mailing list (yade-dev@lists.launchpad.net) as attachment. In either way, after reviewing your changes they will be added to the main trunk.

When the pull request has been reviewed and accepted, your changes are integrated in the main trunk. Everyone will get them via ``git fetch``.

Updating
--------

You may want to get changes done by others::

 git fetch upstream  #Pull new updates from the upstream to your branch. Eq. of "bzr update", updating the remote branch from the upstream yade/trunk
 git merge upstream/master  #Merge upstream changes into your master-branch (eq. of "bzr update", updating your local repository from the remote branch)

Alternatively, this will do fetch+merge all at once (discouraged if you have uncommited changes)::

 git pull

**********************************************************************************************
Working directly on git/trunk (not possible after GitLab migration - section will get updates)
**********************************************************************************************

This direct access to trunk will sound more familiar to `bzr <http://bazaar.canonical.com/en/>`_ or `svn <https://subversion.apache.org/>`_ users. It is only possible for members of the git team "developpers". Send an email at yade-dev@lists.launchpad.net to join this team (don't forget to tell your git account name).

* Get trunk:

   ::

    git clone git@gitlab.com:yade-dev/trunk.git

   This creates a new folder, named trunk, that contains the whole code.

* Update

   ::

    git pull

* Commit to local repository

   ::

    git commit filename1 filename2 ...

* Push changes to remote trunk

   ::

    git push

   Now, the changes you made are included in the on-line code, and can be get back by every user.

   To avoid confusing logs after each commit/pull/push cycle, it is better to setup automatic rebase::

    git config --global branch.autosetuprebase always

   Now your file ~/.gitconfig should include::

	  [branch]
	    autosetuprebase = always

   Check also .git/config file in your local trunk folder (rebase = true)::

	  [remote "origin"]
	    url = git@gitlab.com:yade-dev/trunk.git
	    fetch = +refs/heads/*:refs/remotes/origin/*
	  [branch "master"]
	    remote = origin
	    merge = refs/heads/master
	    rebase = true

   Auto-rebase may have unpleasant side effects by blocking "pull" if you have uncommited changes. In this case you can use "git stash"::

 git pull
 lib/SConscript: needs update
 refusing to pull with rebase: your working tree is not up-to-date
 git stash #hide the uncommited changes away
 git pull  #now it's ok
 git push  #push the commited changes
 git stash pop #get uncommited changes back

********************************************
General guidelines for pushing to yade/trunk
********************************************

1. Set autorebase globaly on the computer (only once see above), or at least on current local branch. Non-rebased pull requests will not be accepted on the upstream. This is to keep history linear, and avoid the merge commits.  

2. Inspect the diff to make sure you will not commit junk code (typically some "cout<<" left here and there), using in terminal:

   ::

    git diff file1
  
   Or using your preferred difftool, such as kdiff3:
  
   ::
  
    git difftool -t kdiff3 file1

   Or, alternatively, any GUI for git: gitg, git-cola... 

3. Commit selectively:

 ::

  git commit file1 file2 file3 -m "message" # is good
  git commit -a -m "message"                # is bad. It is the best way to commit things that should not be commited

4. Be sure to work with an up-to-date version launching:

 ::

  git pull

5. Make sure it compiles and that regression tests pass: try ``yade --test`` and ``yade --check``.

6. You can finally let all Yade-users enjoy your work:

 ::

  git push


**Thanks a lot for your cooperation to Yade!**
