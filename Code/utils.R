
cellCols <- c("iSC"="#e5721c", "pmSC"="#9a5268", "prol. SC"="#ef2700",
              "mSC1"="#4e57b2", "tSC"="#71790d", "mSC2"="#003399",
              "nmSC"="#00af00", "nm(R)SC"="#00af00",
              "mSC3"="#067cfc", "mSC"="#4e57b2")
timeCols <- c("E13.5"="#f38df2", "E17.5"="#4cd3f2",
              "P1"="#d7b73a", "P5"="#7aa1e3",
              "P14"="#65ce8d", "P24"="#f48cc4", "P60"="#adc33c")
timeCols2 <- setNames(timeCols, sub("\\.\\d$", "", names(timeCols)))

bulkCellCols <- c("allSC"="#35918e", "mSC"="#3c3c3b", "nmSC"="#a5a5a5")

cellCols10xP1 <- c("prol. SC"="#66cc66", "iSC"="#079f01", "pmSC"="#006600",
                 "prol. Fb"="#ff9966", "EpC"="#993300", "EnC"="#e84e1b",
                 "PnC"="#ff9900", "FbRel*"="#b2855b", "Per/VSMC"="#1eb1a9",
                 "EC"="#367acc", "IC"="#dd1d20")

cellCols10xP60 <- c("SC"="#079f01","EpC"="#993300","EnC"="#e84e1b",
                    "PnC"="#ff9900","IC"="#e02020",
                    "EC1"="#A57BB4","EC2"="#1439af","Per/EC*"="#3399cc",
                    "Per/VSMC"="#1fb2aa")
