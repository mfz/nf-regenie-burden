# use rclone to deCODE S3 bucket

To be done once in new workspace

- install rclone as user on AllOfUS machine

  #+BEGIN_SRC sh
    mkdir -p ~/bin
    cd /tmp
    curl -O https://downloads.rclone.org/rclone-current-linux-amd64.zip
    unzip rclone-current-linux-amd64.zip
    cd rclone-*-linux-amd64
    cp rclone ~/bin/
    
  #+END_SRC


- transfer rclone.config for ito17499:ito-17499-allofus to AllOfUS machine
  
  ~/.config/rclone/rclone.conf
  
Whenever we need to transfer

- create directory with files to download, and copy files from bucket into it

  #+BEGIN_SRC sh
    DATADIR=mydatadir
    mkdir -p ${DATADIR}
    gsutil -m cp <SRCDIR>/* ${DATADIR}
  #+END_SRC

- create tar archive and split into chunks (AllOfUS does not like files above 100M)

  #+BEGIN_SRC sh
  tar cf ${DATADIR}.tar ${DATADIR}

  mkdir ${DATADIR}chunks
  split -b 90M ${DATADIR}.tar ${DATADIR}chunks/chunk_
  #+END_SRC


- use rclone to transfer chunks to bucket (limit bandwidth so we do not run into egress alarm)
  (Problem: the bandwidth of the AoU workstation is limited by AoU. So the speed is extremely slow)

  #+BEGIN_SRC sh
  rclone copy ./${DATADIR}chunks ito17499:ito-17499-allofus/${DATADIR}chunks --bwlimit 10M
  #+END_SRC


- transfer bucket to deCODE host

  log into office linux

  #+BEGIN_SRC sh
   DATADIR=mydatadir
   rclone copy  ito17499:ito-17499-allofus/${DATADIR}chunks ./${DATADIR}chunks
  #+END_SRC

#+BEGIN_SRC sh
  cat ./${DATADIR}chunks/chunk_* > ${DATADIR}.tar
  tar xf ${DATADIR}.tar
#+END_SRC
