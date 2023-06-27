
# usethis::use_revdep()
# revdepcheck::revdep_check(num_workers = 4)

tt = devtools::spell_check()
View(tt)
devtools::build()
#devtools::release()
