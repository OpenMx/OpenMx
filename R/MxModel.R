

setConstructorS3("MxModel", function() { 

  freeVariablesList <- list();
  
  extend(Object(), "MxModel",
    .freeVariablesList=freeVariablesList
  );

})



setMethodS3("$<-", "MxModel", function(this, name, value) {
  memberAccessorOrder <- attr(this, ".memberAccessorOrder");
  if (is.null(memberAccessorOrder))
    memberAccessorOrder <- c(1,2,3,4,5);

  #
  # This portion of the method is specific to OpenMx
  # Michael Spiegel, May 2, 2007
  #
  if (inherits(value,"MxMatrix")) {
    transformMatrix(value, length(this$.freeVariablesList))
  }

  for (memberAccessor in memberAccessorOrder) {
    if (memberAccessor == 1) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # Search for a set<Name>() method
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      if (is.null(attr(this, "disableSetMethods"))) {
        firstChar <- substr(name, 1,1);
        # Do not try to access private fields using a set<Name>() method,
        # because such a functionality means that the user *expects* that
        # there actually is a field called '.<name>', which he or she
        # should not do since it is a private field!
        # Is it a private field?
        if (!identical(firstChar, ".")) {
          # Field names can not contain spaces...
          if (regexpr(" ", name) == -1) {
            # 1. Is it a set<Name>() method?
            capitalizedName <- name;
            substr(capitalizedName,1,1) <- toupper(firstChar);
            setMethodNames <- paste("set", capitalizedName, ".", class(this), sep="");
            for (setMethodName in setMethodNames) {
              if (exists(setMethodName, mode="function")) {
                ref <- this;
                attr(ref, "disableSetMethods") <- TRUE;
                get(setMethodName, mode="function")(ref, value);
                return(invisible(this));
              }
            }
          } # if ("no space in the name")
        } # if ("is private field")
      } # if (is.null(attr(this, "disableSetMethods")))
    } else if (memberAccessor == 2) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # Search for a <name> field
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # 2. If there exists a field, assign the value to that field.
      envir <- attr(this, ".env");
      if (exists(name, envir=envir)) {
        assign(name, value, envir=envir);

        #
        # This portion of the method is specific to OpenMx
        # Michael Spiegel, May 2, 2007
        #
        if (inherits(value,"MxMatrix")) {
          this$updateFreeVariablesList();
        }
        
        return(invisible(this));
      }
    } else if (memberAccessor == 3) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # Search for a <name> attribute.   /Should this be removed?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # 3. If there exists an attribute field, assign the value to that field.
      if (is.element(name, names(attributes(this)))) {
        attr(this, name) <- value;
        return(invisible(this));
      }
    } else if (memberAccessor == 4) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # Search for a static <name> field
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # 4. If not, it might be that it is a static field
      static <- getStaticInstance(get(class(this)[1]));
      static.envir <- attr(static, ".env");
      # For static method calls, e.g. Object$load, 'this' has no
      # environment assigned and therefore, for now, no static
      # fields.
      if (!is.null(static.envir) && exists(name, envir=static.envir, inherit=FALSE)) {
        assign(name, value, envir=static.envir);
        return(invisible(this));
      }
    } else if (memberAccessor == 5) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # Create a new field <name>
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # 5. Otherwise, assign the value to a new field.
      assign(name, value, envir=envir);
      
      #
      # This portion of the method is specific to OpenMx
      # Michael Spiegel, May 2, 2007
      #
      if (inherits(value,"MxMatrix")) {
        this$updateFreeVariablesList();
      }
            
      return(invisible(this));
    }
  } # for (memberAccessor in memberAccessorOrder)

  invisible(this);
}, createGeneric=FALSE) # $<-()

setMethodS3("[[<-", "MxModel", function(this, name, value) {
  UseMethod("$<-");
#   "$<-"(this, name, value);
}, createGeneric=FALSE) # "[[<-"()

