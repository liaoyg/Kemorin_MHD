/*
//  control_arrays_IO_c.c
//  
//
//  Created by Hiroaki Matsui on 2018/05/24.
*/

#include "control_arrays_IO_c.h"


void alloc_ctl_chara_array(struct chara_ctl_array *c_array){
	int i;
	
	c_array->c_array_item = (struct chara_ctl_item **) malloc(c_array->num*sizeof(struct chara_ctl_item *));
	for(i=0;i<c_array->num;i++){
		c_array->c_array_item[i] = (struct chara_ctl_item *) malloc(sizeof(struct chara_ctl_item));
		alloc_ctl_chara_item(c_array->c_array_item[i]);
	};
	c_array->icou = 0;
	return;
}
void dealloc_ctl_chara_array(struct chara_ctl_array *c_array){
	int i;
	
	if(c_array->num < 0) return;
	for(i=0;i<c_array->num;i++){
		dealloc_ctl_chara_item(c_array->c_array_item[i]);
		free(c_array->c_array_item[i]);
	};
	free(c_array->c_array_item);
	return;
}
void read_character_ctl_array_c(FILE *fp, char *buf, const char *label,
			struct chara_ctl_array *c_array){
	if(c_array->icou > 0) return;
	if(find_control_array_flag_c(buf, label, &c_array->num) == 0) return;
	
	if(c_array->num > 0) alloc_ctl_chara_array(c_array);
	fgets(buf, LENGTHBUF, fp);
	while(find_control_end_array_flag_c(buf, label, c_array->num, c_array->icou) == 0){
		fgets(buf, LENGTHBUF, fp);
		
		read_character_ctl_item_c(buf, label, c_array->c_array_item[c_array->icou]);
		c_array->icou = c_array->icou + c_array->c_array_item[c_array->icou]->iflag;
	};
	
	return;
}
void write_character_ctl_array_c(FILE *fp, int level, int maxlen,
			const char *label, struct chara_ctl_array *c_array){
	int i;
	
	if(c_array->icou == 0) return;
	
	level = write_array_flag_for_ctl_c(fp, level, label, c_array->num);
	for(i=0;i<c_array->num;i++){
		write_character_ctl_item_c(fp, level, strlen(label),
					label, c_array->c_array_item[i]);
	};
	level = write_end_array_flag_for_ctl_c(fp, level, label);
	return;
}


void alloc_ctl_int_array(struct int_ctl_array *i_array){
	int i;
	
	i_array->i_array_item = (struct int_ctl_item **) malloc(i_array->num*sizeof(struct int_ctl_item *));
	for(i=0;i<i_array->num;i++){
		i_array->i_array_item[i] = (struct int_ctl_item *) malloc(sizeof(struct int_ctl_item));
		init_ctl_int_item(i_array->i_array_item[i]);
	};
	i_array->icou = 0;
	return;
}
void dealloc_ctl_int_array(struct int_ctl_array *i_array){
	int i;
	
	if(i_array->num < 0) return;
	for(i=0;i<i_array->num;i++){
		free(i_array->i_array_item[i]);
	};
	free(i_array->i_array_item);
	return;
}
void read_integer_ctl_array_c(FILE *fp, char *buf, const char *label,
			struct int_ctl_array *i_array){
	if(i_array->icou > 0) return;
	if(find_control_array_flag_c(buf, label, &i_array->num) == 0) return;
	
	if(i_array->num > 0) alloc_ctl_int_array(i_array);
	
	fgets(buf, LENGTHBUF, fp);
	while(find_control_end_array_flag_c(buf, label, i_array->num, i_array->icou) == 0){
		if(i_array->icou >= i_array->num){
			printf("Numnber of item is larger than defined \n");
			return;
		}
		
		read_integer_ctl_item_c(buf, label, i_array->i_array_item[i_array->icou]);
		i_array->icou = i_array->icou + i_array->i_array_item[i_array->icou]->iflag;
		
		fgets(buf, LENGTHBUF, fp);
	};
	
	return;
}
void write_integer_ctl_array_c(FILE *fp, int level, int maxlen,
			const char *label, struct int_ctl_array *i_array){
	int i;
	
	if(i_array->icou == 0) return;
	
	level = write_array_flag_for_ctl_c(fp, level, label, i_array->num);
	for(i=0;i<i_array->num;i++){
		write_integer_ctl_item_c(fp, level, strlen(label),
					label, i_array->i_array_item[i]);
	};
	level = write_end_array_flag_for_ctl_c(fp, level, label);
	return;
}


void alloc_ctl_i2_array(struct i2_ctl_array *i2_array){
	int i;
	
	i2_array->i1_array_item = (struct int_ctl_item **) malloc(i2_array->num*sizeof(struct int_ctl_item *));
	i2_array->i2_array_item = (struct int_ctl_item **) malloc(i2_array->num*sizeof(struct int_ctl_item *));
	for(i=0;i<i2_array->num;i++){
		i2_array->i1_array_item[i] = (struct int_ctl_item *) malloc(sizeof(struct int_ctl_item));
		i2_array->i2_array_item[i] = (struct int_ctl_item *) malloc(sizeof(struct int_ctl_item));
		init_ctl_int_item(i2_array->i1_array_item[i]);
		init_ctl_int_item(i2_array->i2_array_item[i]);
	};
	i2_array->icou = 0;
	return;
}
void dealloc_ctl_i2_array(struct i2_ctl_array *i2_array){
	int i;
	
	if(i2_array->num < 0) return;
	for(i=0;i<i2_array->num;i++){
		free(i2_array->i1_array_item[i]);
		free(i2_array->i2_array_item[i]);
	};
	free(i2_array->i1_array_item);
	free(i2_array->i2_array_item);
	return;
}
void read_int2_ctl_array_c(FILE *fp, char *buf, const char *label,
			struct i2_ctl_array *i2_array){
	
	if(i2_array->icou > 0) return;
	if(find_control_array_flag_c(buf, label, &i2_array->num) == 0) return;
	
	if(i2_array->num > 0) alloc_ctl_i2_array(i2_array);
	fgets(buf, LENGTHBUF, fp);
	while(find_control_end_array_flag_c(buf, label, i2_array->num, i2_array->icou) == 0){
		if(i2_array->icou >= i2_array->num){
			printf("Numnber of item is larger than defined \n");
			return;
		}
		
		read_int2_ctl_item_c(buf, label, i2_array->i1_array_item[i2_array->icou],
					i2_array->i2_array_item[i2_array->icou]);
		i2_array->icou = i2_array->icou + i2_array->i1_array_item[i2_array->icou]->iflag;
		
		fgets(buf, LENGTHBUF, fp);
	};
	
	return;
}
void write_int2_ctl_array_c(FILE *fp, int level, int maxlen,
			const char *label, struct i2_ctl_array *i2_array){
	int i;
	
	if(i2_array->icou == 0) return;
	
	level = write_array_flag_for_ctl_c(fp, level, label, i2_array->num);
	for(i=0;i<i2_array->num;i++){
		write_int2_ctl_item_c(fp, level, strlen(label),
					label, i2_array->i1_array_item[i], i2_array->i2_array_item[i]);
	};
	level = write_end_array_flag_for_ctl_c(fp, level, label);
	return;
}
