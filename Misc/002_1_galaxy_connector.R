library(GalaxyConnector)
GX_API_KEY= "0f680c963c66a6a2696231cc8d362f23"    # Your Galaxy API Key
GX_URL= "https://share-galaxy.ibers.aber.ac.uk/"    # The Galaxy instance url


gx_init(GX_API_KEY,GX_URL,gx_list_histories()$id)
gx_list_histories()
gx_list_history_datasets()

gx_current_history()
s <- gx_show_dataset("a7df7729b223b816")
gx_get(4)
