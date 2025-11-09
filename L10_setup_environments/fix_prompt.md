# How to fix the prompt to remove path



## Quick fix in the terminal

```
export PS1='\u@\h$ '
```


## Permanent fix in the codespace


Open (or create) `~/.bashrc` in the Codespace.

Look for a line that sets PS1 (it often includes \w or \W), e.g.:

```
PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
```

Comment it out and add your own minimal prompt, for example:

```
# Old prompt:
# PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '

# New minimal prompt:
PS1='\u@\h$ '
```

Save and reload:

```
source ~/.bashrc
```