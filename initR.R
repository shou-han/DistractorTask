user_name = "shou-han"
user_email = "shn.zhou@gmail.com"

cmd1 <- paste0("git config --global user.name ", user_name)
cmd2 <- paste0("git config --global user.email ", user_email)
cmd3 <- paste0("git config --global --list")
# this is 5 hours
cmd4 <- paste0("git config credential.helper 'cache --timeout=18000'")

system(cmd1)
system(cmd2)
system(cmd3)
system(cmd4)