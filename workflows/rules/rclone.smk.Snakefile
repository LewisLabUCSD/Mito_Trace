rule rclone:
    input:
         "{input}"
    output:
        "{input}/output/rsync_copy.ok"
    params:
        rsync_loc = config["rsync_location"]
    shell: "rcopy {input} drive:{rsync_loc} > rsync_copy.ok && mv rsync_copy.ok {output}/"
