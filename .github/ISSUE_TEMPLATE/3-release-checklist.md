---
name: Release checklist
about: Steps for cutting a new release
title: "Release vX.Y.Z"
---

Release version: vX.Y.Z  <!-- replace everywhere, including the title -->


## Before Release
- [ ] Merge final PRs
   - [ ] #<PR number>
   - [ ] NCAR/ccpp-physics#<PR number>
- [ ] Create list of new feature for release
- [ ] Update tests
- [ ] Update documentation
  - [ ] Update user guides
  - [ ] Check for broken links
- [ ] Update Docker Images
- [ ] Create `Release Notes`

## During Release
- [ ] Create release branch
- [ ] Tag release, push tag if created with CLI
- [ ] Draft GitHub release
- [ ] Create tarball of code with submodules (Github doesn't include them in their tarball)
```bash
$ git clone --branch vX.Y.Z --depth 1 \
        --recurse-submodules --shallow-submodules \
        https://github.com/NCAR/ccpp-scm.git ccpp-scm-vX.Y.Z
$ git -C ccpp-scm-vX.Y.Z submodule status --recursive
$ tar -czf vX.Y.Z.tar.gz \
    --exclude='.git' \
    --exclude='.DS_Store'  --exclude='._*'  --exclude='__MACOSX' \
    --owner=0 --group=0 --numeric-owner \
    ccpp-scm-vX.Y.Z
$ sha256sum vX.Y.Z.tar.gz > vX.Y.Z.tar.gz.sha256
```
- [ ] Create zip file with submodules, after following instructions for tarball
```bash
$ zip -r -X -y vX.Y.Z.zip ccpp-scm-vX.Y.Z \
    -x '*/.git' '*/.git/*' \
       '*/.DS_Store' '*/._*' \
       '*/__MACOSX' '*/__MACOSX/*'
$ sha256sum vX.Y.Z.zip > vX.Y.Z.zip.sha256
```
- [ ] Add artifacts to release
   - [ ] Add tarball
   - [ ] Add zip file
   - [ ] See [v8.0.0 release](https://github.com/NCAR/ccpp-scm/releases#release-v8.0.0) for full list
   - [ ] tutorial_files.tar.gz
- [ ] [DTC website](https://dtcenter.org/software-tools/common-community-physics-package-ccpp/download): update information from previous `vX.Y.Z` to new `vX.Y.Z`
- [ ] Tag Framework [ccpp-scm_vX.Y.Z](https://github.com/NCAR/ccpp-framework/releases/tag/ccpp-scm_vX.Y.Z)
- [ ] Tag Techdoc [vX.Y.Z](https://github.com/NCAR/ccpp-doc) and add [readthedocs version](https://app.readthedocs.org/projects/ccpp-techdoc/)

## After release
- [ ] Verify new Zenodo DOI
- [ ] Update training lessons to `vX.Y.Z`
- [ ] Make sure CCPP-Physics release also occurred
- [ ] Push Docker image `vX.Y.Z`
