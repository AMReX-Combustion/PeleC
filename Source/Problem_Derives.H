// problem-specific derived variables would be put here

//
// Error between solution and exact solution for MMS
//
if (do_mms){
  derive_lst.add("rhommserror",IndexType::TheCellType(),1,pc_derrhommserror,the_same_box);
  derive_lst.addComponent("rhommserror",desc_lst,State_Type,Density,NUM_STATE);

  derive_lst.add("ummserror",IndexType::TheCellType(),1,pc_derummserror,the_same_box);
  derive_lst.addComponent("ummserror",desc_lst,State_Type,Density,NUM_STATE);

  derive_lst.add("vmmserror",IndexType::TheCellType(),1,pc_dervmmserror,the_same_box);
  derive_lst.addComponent("vmmserror",desc_lst,State_Type,Density,NUM_STATE);

  derive_lst.add("wmmserror",IndexType::TheCellType(),1,pc_derwmmserror,the_same_box);
  derive_lst.addComponent("wmmserror",desc_lst,State_Type,Density,NUM_STATE);

  derive_lst.add("pmmserror",IndexType::TheCellType(),1,pc_derpmmserror,the_same_box);
  derive_lst.addComponent("pmmserror",desc_lst,State_Type,Density,NUM_STATE);
 }
