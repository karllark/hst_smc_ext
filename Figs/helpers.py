

def prettyname(name):

    pname = None
    if "2dfs" in name:
        pname = name.replace("2dfs", "2dFS ")
        pname = pname.replace(" 0", " ")
    if "azv" in name:
        pname = name.replace("azv", "AzV ")
    if "bbb" in name:
        pname = name.replace("bbb-smc", "BBB SMC")
    if "mr12" in name:
        pname = name.replace("mr12-star", "MR12 ")
    if "ngc330" in name:
        pname = name.replace("ngc330-", "NGC330 ELS ")
    if "ngc346" in name:
        pname = name.replace("ngc346-", "NGC330 ELS ")
    if "sk" in name:
        pname = name.replace("sk", "SK ")
    if "smc5" in name:
        pname = name.replace("smc5-", "SMC5-")
        pname = pname.replace("0", "")
    if "smcave" in name:
        pname = "SMC Average"
    if "smcave_bumps" in name:
        pname = "SMC Bumps"

    return pname